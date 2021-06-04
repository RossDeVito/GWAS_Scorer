import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.preprocessing import StandardScaler


if __name__ == '__main__':
    log_score = True

    # Load Fst data
    Fst_df = pd.read_csv('data/preprocessed/Fst_v1.csv')
    Fst_dict = Fst_df.set_index(['pop_1', 'pop_2']).to_dict()

    Fst_dict['hudson_fst'][('CEU', 'YRI')]

    # Load study scores
    save_dir = 'data/results'

    height_studies = [
        'GCST005908',   #EUR1
        'GCST000817',   #EUR2
        'GCST000611',   #JAP1
        'GCST008839',   #JAP2
        'GCST001263',   #AFR1
        'GCST001290',	#AFR2
        'GCST002702',   #EAS1
        'GCST004212',   #GBR1
        'GCST008163',   #EUR3
        'GCST002647',	#EUR4
        'GCST001885',   #CHB1
        'GCST001956',   #EUR5
        'GCST000372',	#EUR6
        'GCST000174',	#EUR7
        'GCST000175',	#EUR8
        'GCST000176',	#EUR9
    ]

    schizophrenia_studies = [
        'GCST006803',   #EUR1_S
        'GCST003880',   #CHB1_S
        'GCST001242',	#EUR2_S
        'GCST007205',	#JPT1_S
        'GCST001757',	#EUR3_S
        'GCST002149',	#EUR4_S
        'GCST001657',	#EUR5_S
    ]

    bmi_studies = [
        'GCST009004',	#EUR1_B
        'GCST009003',	#EUR1B_B
        'GCST009001',	#EUR1c_B
    ]

    panc_cancer_studies = [
        'GCST005786',	#AFR1_PC
    ]

    studies = panc_cancer_studies

    study_to_sample_pop = {
        'GCST005908': 'EUR',    #EUR1
        'GCST000817': 'EUR',    #EUR2
        'GCST000611': 'JPT',    #JAP1
        'GCST008839': 'JPT',    #JAP2
        'GCST001263': 'AFR',    #AFR1
        'GCST001290': 'AFR',	#AFR2
        'GCST002702': 'EAS',    #EAS1
        'GCST004212': 'GBR',	#GBR1
        'GCST008163': 'EUR',    #EUR3
        'GCST002647': 'EUR',    #EUR4
        'GCST001885': 'CHB',    #CHB1
        'GCST001956': 'EUR',    #EUR5
        'GCST000372': 'EUR',	#EUR6
        'GCST000174': 'EUR',	#EUR7
        'GCST000175': 'EUR',	#EUR8
        'GCST000176': 'EUR',	#EUR9
        'GCST006803': 'EUR',    #EUR1_S
        'GCST003880': 'CHB',    #CHB1_S
        'GCST001242': 'EUR',	#EUR2_S
        'GCST007205': 'JPT',	#JPT1_S
        'GCST001757': 'EUR',	#EUR3_S
        'GCST002149': 'EUR',	#EUR4_S
        'GCST001657': 'EUR',	#EUR5_S
        'GCST009004': 'EUR',	#EUR1_B
        'GCST009001': 'EUR',	#EUR1c_B
        'GCST009003': 'EUR',	#EUR1B_B
        'GCST005786': 'AFR',	#AFR1_PC
    }

    # add sample pops Fst with self as 0
    for pop in np.unique(list(study_to_sample_pop.values())):
        Fst_dict['hudson_fst'][(pop, pop)] = 0
        Fst_dict['patterson_fst'][(pop, pop)] = 0

    # load data and add scaled score
    all_dfs = []
    std_scaler = StandardScaler()

    for s in studies:
        df = pd.read_csv(os.path.join(save_dir, '{}_scores.csv'.format(s)))

        if log_score:
            df['scores'] = np.log(df.scores)

        df['scaled_score'] = std_scaler.fit_transform(df.scores.values.reshape(-1,1))
        df['gender_scaled_score'] = 0
        df.loc[df.gender == 'male', 'gender_scaled_score'] = std_scaler.fit_transform(
            df.scores[df.gender == 'male'].values.reshape(-1,1)
        )
        df.loc[df.gender == 'female', 'gender_scaled_score'] = std_scaler.fit_transform(
            df.scores[df.gender == 'female'].values.reshape(-1,1)
        )

        df['sample_pop'] = study_to_sample_pop[s]
        df['study'] = df[['STUDY ACCESSION', 'sample_pop']].agg(' '.join, axis=1)

        all_dfs.append(df)

    df = pd.concat(all_dfs)

    # Remove samples with null pop labels
    df = df.dropna(subset=['pop', 'super_pop'])

    # Add Fst
    df['pop_Fst'] = df.apply(lambda s: Fst_dict['hudson_fst'][(s['pop'], s['sample_pop'])], axis=1)
    df['super_pop_Fst'] = df.apply(lambda s: Fst_dict['hudson_fst'][(s['super_pop'], s['sample_pop'])], axis=1)

    # Add in study variance by raw and std score
    df['pop_score_variance'] = df.groupby(['STUDY ACCESSION', 'pop']).scores.transform('var')
    df['pop_score_variance_bg'] = df.groupby(['STUDY ACCESSION', 'pop', 'gender']).scores.transform('var')

    # df['super_pop_score_variance'] = df.groupby(['STUDY ACCESSION', 'super_pop', 'gender']).scores.transform('var')
    # df['super_pop_scaled_score_variance'] = df.groupby(['STUDY ACCESSION', 'super_pop', 'gender']).scaled_score.transform('var')
    # df['super_pop_gender_scaled_score_variance'] = df.groupby(['STUDY ACCESSION', 'super_pop', 'gender']).gender_scaled_score.transform('var')

    # df['pop_score_median'] = df.groupby(['STUDY ACCESSION', 'pop', 'gender']).scores.transform('median')
    # df['pop_scaled_score_median'] = df.groupby(['STUDY ACCESSION', 'pop', 'gender']).scaled_score.transform('median')
    # df['pop_gender_scaled_score_median'] = df.groupby(['STUDY ACCESSION', 'pop', 'gender']).gender_scaled_score.transform('median')

    # Plot
    with sns.plotting_context(context='paper', font_scale=1.):
        # Plot variances
        # sns.relplot(data=df, x='pop_Fst', y='pop_gender_scaled_score_variance',
        #             hue='super_pop', style='pop', col='gender',
        #             estimator=np.median,
        #             facet_kws={'margin_titles':True})
        # plt.tight_layout()
        # plt.show()

        # sns.relplot(data=df, x='pop_Fst', y='pop_gender_scaled_score_variance',
        #             hue='super_pop', style='pop',
        #             row='gender', col='study',
        #             estimator=np.median,
        #             facet_kws={'margin_titles':True})
        # plt.tight_layout()
        # plt.show()

        # sns.relplot(data=df, x='pop_Fst', y='pop_gender_scaled_score_variance',
        #             hue='super_pop', style='pop',
        #             row='gender', col='sample_pop',
        #             estimator=np.median,
        #             facet_kws={'margin_titles':True})
        # plt.tight_layout()
        # plt.show()

        sns.relplot(data=df, x='pop_Fst', y='pop_score_variance_bg',
                    hue='super_pop', style='pop',
                    row='gender', col='study',
                    estimator=np.median,
                    facet_kws={'margin_titles':True,
                                'sharey':False})
        plt.tight_layout()
        plt.show()

        # sns.relplot(data=df, x='pop_Fst', y='pop_scaled_score_variance',
        #             hue='super_pop', style='pop',
        #             col='study', estimator=np.median,
        #             facet_kws={'margin_titles':True})
        # plt.tight_layout()
        # plt.show()

        sns.relplot(data=df, x='pop_Fst', y='pop_score_variance',
                    hue='super_pop', style='pop',
                    col='study', estimator=np.median,
                    facet_kws={'margin_titles':True,
                                'sharey':False})
        plt.tight_layout()
        plt.show()

        ## Plot medians
        # sns.relplot(data=df, x='pop_Fst', y='pop_gender_scaled_score_median',
        #             hue='super_pop', style='pop', col='gender',
        #             facet_kws={'margin_titles':True})
        # plt.tight_layout()
        # plt.show()

        # sns.relplot(data=df, x='pop_Fst', y='pop_gender_scaled_score_median',
        #             hue='super_pop', style='pop',
        #             row='gender', col='study',
        #             facet_kws={'margin_titles':True})
        # plt.tight_layout()
        # plt.show()

        # sns.relplot(data=df, x='pop_Fst', y='pop_gender_scaled_score_median',
        #             hue='super_pop', style='pop',
        #             row='gender', col='sample_pop',
        #             facet_kws={'margin_titles':True})
        # plt.tight_layout()
        # plt.show()

        # sns.relplot(data=df, x='pop_Fst', y='pop_score_median',
        #             hue='super_pop', style='pop',
        #             row='gender', col='study',
        #             facet_kws={'margin_titles':True,
        #                         'sharey':False})
        # plt.tight_layout()
        # plt.show()

        # sns.relplot(data=df, x='pop_Fst', y='pop_score_median',
        #             hue='super_pop', style='pop',
        #             row='gender', col='sample_pop',
        #             facet_kws={'margin_titles':True,
        #                         'sharey':False})
        # plt.tight_layout()
        # plt.show()

        # sns.relplot(data=df, x='pop_Fst', y='pop_scaled_score_median',
        #             hue='super_pop', style='pop',
        #             col='study',
        #             facet_kws={'margin_titles':True})
        # plt.tight_layout()
        # plt.show()

        # sns.relplot(data=df, x='pop_Fst', y='pop_score_median',
        #             hue='super_pop', style='pop',
        #             col='study',
        #             facet_kws={'margin_titles':True,
        #                         'sharey':False})
        # plt.tight_layout()
        # plt.show()

