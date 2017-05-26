import cobra, scipy, argparse, pandas
from scipy import Inf

pandas.set_option('display.max_colwidth', 9999999)  # otherwise very long ids are truncated


model = cobra.io.read_sbml_model('/home/cossio/W/2017/DFBAeq/models/mus-musculus/hefzi2016/Supp_Data_S1-iCHOv1/iCHOv1_K1_final.xml')
rxns = pandas.DataFrame([[rxn.id, rxn.lower_bound, rxn.upper_bound] for rxn in model.reactions],
                        columns = ['id', 'lb', 'ub'])
mets = pandas.DataFrame({'id': [met.id for met in model.metabolites]})

# biomass, mmol/gDW
mets['y'] = 0.
for met in model.reactions.biomass_cho.metabolites:
    mets.loc[mets['id'] == met.id, 'y'] = -model.reactions.biomass_cho.get_coefficient(met)

# maintenance demand. Just an ATP hydrolysis.
# Kilburn1969 reports 17 pmol/cell/day. Use CHO dry weight (315 pg, Hefzi2016) to obtain the value used here (in mmol/gDW/h)
mets['e'] = 0.
for met in model.reactions.DM_atp_c_.metabolites:
    mets.loc[mets['id'] == met.id, 'e'] = -model.reactions.DM_atp_c_.get_coefficient(met) * 2.24868
model.reactions.DM_atp_c_.lower_bound = model.reactions.DM_atp_c_.upper_bound = 0

# flux costs
C = 0.078  # Shlomi2011. Solvent capacity. Unitless, since this is the mass fraction occupied by enzymes.
Shlomi = pandas.read_table('Shlomi.txt', dtype={'id': str, 'kcat_fwd': float, 'kcat_bkwd': float, 'mw_fwd': float, 'mw_bkwd': float})
# missing values are replaced by the median of known values
ShlomiMissingKcat = 25. * 3600.   # factor 3600 converts 1/s to 1/hr
ShlomiMissingMW = 85300. / 1000.  # factor 1/1000 converts g/mol to g/mmol

rxns['ap'] = rxns['an'] = ShlomiMissingMW / ShlomiMissingKcat / C

for rxn in model.reactions:
    row = Shlomi[Shlomi['id'] == rxn.id]
    if not row.empty:
        kcat_fwd = row['kcat_fwd'].iloc[0] * 3600.
        kcat_bkwd = row['kcat_bkwd'].iloc[0] * 3600.
        mw_fwd = row['mw_fwd'].iloc[0] / 1000.
        mw_bkwd = row['mw_bkwd'].iloc[0] / 1000.

        if scipy.isfinite(kcat_fwd) and scipy.isfinite(mw_fwd):
            rxns.loc[rxns['id'] == rxn.id, 'ap'] = mw_fwd / kcat_fwd / C
        if scipy.isfinite(kcat_bkwd) and scipy.isfinite(mw_bkwd):
            rxns.loc[rxns['id'] == rxn.id, 'an'] = mw_bkwd / kcat_bkwd / C

assert all(rxns['ap'] >= 0) and all(rxns['an'] >= 0)

mets['bp'] = mets['bn'] = 0.

mets['L'] = mets['V'] = 0.
for rxn in model.reactions:
    if rxn.id.startswith('EX_'):
        for met in rxn.reactants:
            mets.loc[mets['id'] == met.id, 'V'] = Inf
            if rxn.upper_bound > 0 or rxn.upper_bound == rxn.lower_bound:
                mets.loc[mets['id'] == met.id, 'L'] = -Inf
            else:
                mets.loc[mets['id'] == met.id, 'V'] = 0.

# limit uptake of substrates
mets.loc[mets['id'] == 'glc_D_e', 'V'] = Vglc = 0.5  # mmol/gDW/h, Nolan2011, Kiparissides2011
for met in {'arg_L_e', 'cys_L_e', 'gln_L_e', 'his_L_e', 'ile_L_e', 'leu_L_e', 'lys_L_e',
            'met_L_e', 'phe_L_e', 'thr_L_e', 'trp_L_e', 'tyr_L_e', 'val_L_e',
            'gthrd_e', '4hpro_LT_e', 'btn_e', 'fol_e', 'inost_e', 'ncam_e', 'ribflv_e'}:
    mets.loc[mets['id'] == met, 'V'] = Vglc / 10.

# set to zero reactions that consume/exchange stuff
for rxn in model.reactions:
    if not rxn.products or not rxn.reactants:
        if not rxn.id.startswith('SK_'):
            # SK_ rxns are needed by Hefzi2016 model to live
            rxns.loc[rxns['id'] == rxn.id, 'lb'] = rxns.loc[rxns['id'] == rxn.id, 'ub'] = 0.

sto = pandas.DataFrame([[pos[0] + 1, pos[1] + 1, v]
                        for (pos, v) in model.to_array_based_model().S.asformat(format='dok').iteritems()],
                       columns=['i', 'k', 's'])
sto['met'] = [model.metabolites[i - 1].id for i in sto['i']]
sto['rxn'] = [model.reactions[k - 1].id   for k in sto['k']]


out = 'cho'

with open(out + '.sto', 'w') as f:
    f.write(sto.to_string(index=False))
with open(out + '.rxns', 'w') as f:
    f.write(rxns.to_string(index=False))
with open(out + '.mets', 'w') as f:
    f.write(mets.to_string(index=False))


print 'm =', len(model.metabolites), 'n =', len(model.reactions)
