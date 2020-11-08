# Functions based on code from LFPtools by jsosa:
import pandas as pd

def write_bdy(bdylfp, runcsv, t):
    """
    Subroutine used to write .BDY file
    Inflows are based on JRC hydrological model output
    """

    print("     writing .bdy file...",bdylfp)

    run = pd.read_csv(runcsv, index_col=0)

    # Select only date columns
    rund = run[[i for i in run.columns if (i[0] == '1') | (i[0] == '2')]].T

    # creating file
    with open(bdylfp, 'w') as f:
        f.write('# GBM-flood bdy file'+'\n')

    # writing inflows
    for i in rund.columns:
        try:
            r = rund[i].to_frame()
        except:
            print('Error, this is probably because multiple rows have the same linkID',rund[i])
            raise
        #    r = rund[i].copy()
        #print(i,type(r))
        #r = r.to_frame()
        #print(type(r))
        r['hours'] = list(range(0, t*24, 24))
        with open(bdylfp, 'a') as f:
            f.write('in'+str(i)+'\n')
            f.write(str(len(r['hours']))+' '+'hours'+'\n')
        r.to_csv(bdylfp, sep=' ', float_format='%.7f',
                 index=False, header=False, mode='a')


def write_bci(bcilfp, runcsv):
    """
    Writes bcif: XXX.bci file to be used in LISFLOOD-FP
    Uses runfcsv: XXX_run.csv
    """

    print("     writing .bci file...")

    run = pd.read_csv(runcsv, index_col=0)

    runi = run[['x', 'y']].T

    # creating file
    with open(bcilfp, 'w') as f:
        f.write('# euflood bci file'+'\n')

    # writing inflows
    with open(bcilfp, 'a') as f:
        for i in runi.columns:
            t = 'P'
            x = str(runi[i].loc['x'])
            y = str(runi[i].loc['y'])
            n = 'in' + str(i)
            f.write(t + ' ' + x + ' ' + y + ' ' + 'QVAR' + ' ' + n + '\n')
    
    # Writing other bundary conditions (eg. wall)

    with open(bcilfp,'a') as f:
        f.write('N -9999 9999 FREE' + '\n')
        f.write('S -9999 9999 FREE' + '\n')
        f.write('E -9999 9999 FREE' + '\n')
        f.write('W -9999 9999 FREE' + '\n')

def write_bci_v2(bcilfp, runcsv,downbc):
    """
    Writes bcif: XXX.bci file to be used in LISFLOOD-FP
    Uses runfcsv: XXX_run.csv
    """

    print("     writing .bci file...")

    run = pd.read_csv(runcsv, index_col=0)

    runi = run[['x', 'y']].T

    # creating file
    with open(bcilfp, 'w') as f:
        f.write('# euflood bci file'+'\n')

    # writing inflows
    with open(bcilfp, 'a') as f:
        for i in runi.columns:
            t = 'P'
            x = str(runi[i].loc['x'])
            y = str(runi[i].loc['y'])
            n = 'in' + str(i)
            f.write(t + ' ' + x + ' ' + y + ' ' + 'QVAR' + ' ' + n + '\n')
    
    # Writing other bundary conditions (eg. wall)
        for bc in downbc:
            f.write(bc+ '\n')
