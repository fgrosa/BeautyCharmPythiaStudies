import argparse
import math
import numpy as np
import uproot
import pandas as pd
import matplotlib.pyplot as plt
import yaml


# Definition of label and colors
Blabels = {
    511: r'B$^{0}$',
    521: r'B$^{+}$',
    531: r'B$_{\rm s}^{0}$',
    5122: r'L$_{\rm b}^{+}$'
}
Bnames = {
    511: 'Bzero',
    521: 'Bplus',
    531: 'Bszero',
    5122: 'Lb',
}
Dlabels = {
    421: r'D$^{0}$',
    -421: r'$\overline{\rm D}^{0}$',
    431: r'D$_{\rm s}^{+}$',
    -431: r'D$_{\rm s}^{-}$',
    411: r'D$^{+}$',
    -411: r'D$^{-}$'
}
Dcol = {
    421: 'red',
    -421: 'darkred',
    431: 'darkorange',
    -431: 'gold',
    411: 'forestgreen',
    -411: 'teal'
}
Dnames = {
    421: 'Dzero',
    -421: 'Dzerobar',
    431: 'Dsplus',
    -431: 'Dsminus',
    411: 'Dplus',
    -411: 'Dminus'
}


def flatten(df, columns):
    idx = df.index.repeat(df[columns[0]].str.len())
    dfflat = pd.concat([
        pd.DataFrame({x: np.concatenate(df[x].values)}) for x in columns], axis=1)
    dfflat.index = idx

    return dfflat.join(df.drop(columns, 1), how='left')


def main():
    '''
    Main function
    '''
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument('configfile', metavar='text', default='config_nonpromptD_diffBR.yml',
                        help='input yaml config file name')
    parser.add_argument('inputfile', metavar='text', default='Bhadrons_diffBRtoD.root',
                        help='input root file name')
    args = parser.parse_args()

    with open(args.configfile, 'r') as ymlCfgFile:
        cfg = yaml.load(ymlCfgFile, yaml.FullLoader)
    Bspecies = cfg['pdgCodeB']
    Dspecies = cfg['pdgCodeD']

    Blabelstouse = {Bhad:Blabels[Bhad] for Bhad in Bspecies}
    Dlabelstouse = {Dmes:Dlabels[Dmes] for Dmes in Dspecies}
    Dcoltouse = {Dmes:Dcol[Dmes] for Dmes in Dspecies}

    tree = uproot.open(args.inputfile)['fTreeDecays']
    df = tree.pandas.df(flatten=False)

    dfB, dfBsel = {}, {}

    for Bhad in Bspecies:
        dfB[Bhad] = df.query('pdgB == {0}'.format(Bhad)).loc[:, ['pD', 'yD', 'pdgD']]
        dfBsel[Bhad] = {}
        for Dmes in Dspecies:
            dfBsel[Bhad][Dmes] = dfB[Bhad][dfB[Bhad]['pdgD'].apply(lambda pdgs: Dmes in pdgs)]
            dfBsel[Bhad][Dmes] = flatten(dfBsel[Bhad][Dmes], ['pD', 'yD', 'pdgD']).query('pdgD == {0}'.format(Dmes))

    BRint, BRintunc, BR, BRunc, BRdf, fig, figint = ({} for iDic in range(7))

    deltap = 0.01
    pmins = [ipt*deltap for ipt in range(0, int(3/deltap))]
    pmaxs = [ipt*deltap for ipt in range(1, (int(3/deltap)+1))]
    plims = list(pmins)
    plims.append(pmaxs[-1])

    # integrated BR
    for Bhad in Bspecies:
        BRint[Bhad] = {}
        BRintunc[Bhad] = {}
        figint[Bhad] = plt.figure(figsize=[7, 7])
        for Dmes in Dspecies:
            br = len(dfBsel[Bhad][Dmes]) / len(dfB[Bhad])
            BRint[Bhad][Dmes] = br
            BRintunc[Bhad][Dmes] = (br * (1-br)) / len(dfB[Bhad])

        plt.grid()
        plt.scatter(list(Dlabelstouse.values()), list(
            BRint[Bhad].values()), color=list(Dcoltouse.values()))
        plt.ylim((1.e-3, 1.))
        plt.yscale('log')
        plt.ylabel(
            r'BR({0} $\rightarrow$ D+X)'.format(Blabelstouse[Bhad]), fontsize=14)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.tight_layout()
        figint[Bhad].savefig('BRint{0}.pdf'.format(Bnames[Bhad]))

    # p* differential BR
    for Bhad in Bspecies:
        BR[Bhad] = {}
        BRunc[Bhad] = {}
        BRdf[Bhad] = {}
        fig[Bhad] = plt.figure(figsize=[7, 7])
        for Dmes in Dspecies:
            BR[Bhad][Dmes], _ = np.histogram(dfBsel[Bhad][Dmes]['pD'], bins=plims)
            BRunc[Bhad][Dmes] = np.sqrt(BR[Bhad][Dmes])
            BR[Bhad][Dmes] = np.divide(BR[Bhad][Dmes], len(dfB[Bhad]))
            BRunc[Bhad][Dmes] = np.divide(BRunc[Bhad][Dmes], len(dfB[Bhad]))

            BRdf[Bhad][Dmes] = pd.DataFrame(list(zip(pmins, pmaxs, BR[Bhad][Dmes], \
                BRunc[Bhad][Dmes])), columns=['pstar_min', 'pstar_max', 'BR', 'BRunc'])
            BRdf[Bhad][Dmes].to_csv(r'{0}to{1}.txt'.format(
                Bnames[Bhad], Dnames[Dmes]), index=None, sep=' ', float_format='%.8f')

            pcent = (np.array(plims)[1:] + np.array(plims)[:-1]) / 2

            plt.subplots_adjust(left=0.13, right=0.98, top=0.95)
            plt.hist(pcent, bins=len(BR[Bhad][Dmes]), weights=np.array(BR[Bhad][Dmes]),
                     range=(min(plims), max(plims)), alpha=0.25, color=Dcoltouse[Dmes],
                     label=r'{0} $\rightarrow$ {1} + X'.format(
                         Blabelstouse[Bhad], Dlabelstouse[Dmes]),
                     histtype='stepfilled', log=True)
            plt.hist(pcent, bins=len(BR[Bhad][Dmes]), weights=np.array(BR[Bhad][Dmes]),
                     range=(min(plims), max(plims)), color=Dcoltouse[Dmes], histtype='step', log=True)
            plt.xlabel('p* (GeV/c)', fontsize=14)
            plt.ylabel(r'BR({0} $\rightarrow$ D+X) per {1} GeV/$c$'.format(Blabelstouse[Bhad], \
                deltap), fontsize=14)
            plt.ylim([1.e-5, 1.])
            plt.xticks(fontsize=14)
            plt.yticks(fontsize=14)
            plt.legend(loc='best')
            fig[Bhad].savefig('BRvspstar{0}.pdf'.format(Bnames[Bhad]))

    plt.show()


main()
