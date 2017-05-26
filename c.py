import cobra, scipy, argparse, pandas
from scipy import Inf, NaN

pandas.set_option('display.max_colwidth', 9999999)  # otherwise very long ids are truncated

model = cobra.io.read_sbml_model('/home/cossio/W/2017/DFBAeq/models/mus-musculus/hefzi2016/Supp_Data_S1-iCHOv1/iCHOv1_K1_final.xml')

# RPMI-1640, Moore1967
IMDM = pandas.read_table('/home/cossio/W/2017/DFBAeq/models/mus-musculus/hefzi2016/CHO-K1/media/IMDM.txt',
                         comment='#', sep='\s+')
media = pandas.DataFrame([[met.id, 0.] for met in model.metabolites], columns=['id', 'mM'])

for met in IMDM['id']:
    assert met in model.metabolites
    q = IMDM[IMDM['id'] == met]
    if not q.empty:
        assert len(q) == 1
        media.loc[media['id'] == met, 'mM'] = q['mM'].iloc[0]

# stuff with infinite supply
for met in {'h2o_e', 'h_e', 'o2_e', 'na1_e', 'k_e', 'ca2_e', 'cl_e', 'so4_e', 'fe2_e', 'fe3_e'}:
    assert met in model.metabolites
    media.loc[media['id'] == met, 'mM'] = Inf

with open('c.txt', 'w') as f:
    f.write(media.to_string(index=False))
