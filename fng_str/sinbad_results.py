import numpy as np
import pandas as pd

'''
parser that reads and turns into pandas dataframes experimental and mcnp results provided 
with sinbad
'''

def read_measured_data(data_filename, sheet_name):
    df = pd.read_excel(data_filename, sheet_name=sheet_name)
    df[['y (cm)', 'E', 'Error on E']] = df['y (cm)      E      Error on E'].str.split('     ', 2, expand=True)
    df = df.iloc[:, 1:]
    try:
        df[:] = df[:].astype(float)
        df['Error on E'] = df['Error on E']*df['E']/100
    except (ValueError, TypeError):
        pass

    return df

def read_measured_data_wpos(data_filename, sheet_name):
    df = pd.read_excel(data_filename, sheet_name=sheet_name)
    df[['Position/y (cm)', 'E', 'Error on E']] = df['Position/y (cm)     E     Error on E'].str.split('    ', 2, expand=True)
    df = df.iloc[:, 1:]
    for column in df:
        try:
            df[column] = df[column].astype(float)
        except (ValueError, TypeError):
            pass
    
    try:
        df['Error on E'] = df['Error on E']*df['E']/100
    except (ValueError, TypeError):
        pass
    
    return df

def read_computed_data(data_filename, sheet_name):

    df = pd.read_excel(data_filename, sheet_name=sheet_name)

    try:
        df[['y(cm)', 'C EFF-3', 'EFF-3 err', 'C FENDL-1', 'FENDL-1 err', 'C FENDL-2', 'FENDL-2 err', 'C/E EFF-3', 'C/E FEN-1', 'C/E FEN-2']] = df['y(cm)     C EFF-3     EFF-3 err     C FENDL-1     FENDL-1 err     C FENDL-2     FENDL-2 err     C/E EFF-3     C/E FEN-1     C/E FEN-2'].str.split('     ', 10, expand=True)
    except KeyError:
        df[['y(cm)', 'C EFF-3', 'EFF-3 err', 'C FENDL-1', 'FENDL-1 err', 'C/E EFF-3', 'C/E FEN-1']] = df['y(cm)     C EFF-3     EFF-3 err     C FENDL-1     FENDL-1 err     C/E EFF-3     C/E FEN-1'].str.split('     ', 7, expand=True)
    
    df = df.iloc[:, 1:]

    for column in df:
        try:
            df[column] = df[column].astype(float)
        except (ValueError, TypeError):
            pass

    try:
        df['EFF-3 err'] = df['EFF-3 err']*df['C EFF-3']/100
        df['FENDL-1 err'] = df['FENDL-1 err']*df['C FENDL-1']/100
        df['FENDL-2 err'] = df['FENDL-2 err']*df['C FENDL-2']/100
    except (ValueError, TypeError, KeyError):
        pass

    return df

def read_computed_data_wpos(data_filename, sheet_name):

    df = pd.read_excel(data_filename, sheet_name=sheet_name)

    try:
        df[['Position/y(cm)', 'C EFF-3', 'EFF-3 err', 'C FENDL-1', 'FENDL-1 err', 'C FENDL-2', 'FENDL-2 err', 'C/E EFF-3', 'C/E FEN-1', 'C/E FEN-2']] = df['Position/y(cm)     C EFF-3     EFF-3 err     C FENDL-1     FENDL-1 err     C FENDL-2     FENDL-2 err     C/E EFF-3     C/E FEN-1     C/E FEN-2'].str.split('     ', 10, expand=True)
    except KeyError:
        df[['Position/y(cm)', 'C EFF-3', 'EFF-3 err', 'C FENDL-1', 'FENDL-1 err', 'C/E EFF-3', 'C/E FEN-1']] = df['Position/y(cm)     C EFF-3     EFF-3 err     C FENDL-1     FENDL-1 err     C/E EFF-3     C/E FEN-1'].str.split('     ', 7, expand=True)
    
    df = df.iloc[:, 1:]

    for column in df:
        try:
            df[column] = df[column].astype(float)
        except (ValueError, TypeError):
            pass
    
    try:
        df['EFF-3 err'] = df['EFF-3 err']*df['C EFF-3']/100
        df['FENDL-1 err'] = df['FENDL-1 err']*df['C FENDL-1']/100
        df['FENDL-2 err'] = df['FENDL-2 err']*df['C FENDL-2']/100
    except (ValueError, TypeError, KeyError):
        pass

    return df