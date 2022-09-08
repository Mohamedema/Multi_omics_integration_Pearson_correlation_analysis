import pandas as pd
main= pd.read_csv("/Users/emam22/Desktop/MOFA/required_file_DL.txt", sep='\t', index_col ="Class")
main_control=main.loc['Control']
main_pre=main.loc['Prediabetic']
main_pre.to_csv("/Users/emam22/Desktop/prediabetic.csv",index=False)
main_control.to_csv("/Users/emam22/Desktop/control.csv",index=False)
#---------------------------------------------------
#Control_view
con= pd.read_csv("/Users/emam22/Desktop/control.csv", sep=',')
#Deseq
first_view_con=con.iloc[:,0:784]
control_sample_id=con.iloc[:,:1]
#metabolome
second_view_con=con.iloc[:,784:918]
second_view_con= control_sample_id.join(second_view_con)
#proteome
third_view_con=con.iloc[:,918:976]
third_view_con=control_sample_id.join(third_view_con)
#gut
forth_view_con=con.iloc[:,976:983]
forth_view_con=control_sample_id.join(forth_view_con)
#nares
fifth_view_con=con.iloc[:,983:]
fifth_view_con=control_sample_id.join(fifth_view_con)
#---------------------------------------------------
#Prediabetic_view
pre= pd.read_csv("/Users/emam22/Desktop/prediabetic.csv", sep=',')
#Deseq
first_view_pre=pre.iloc[:,0:784]
pre_sample_id=pre.iloc[:,:1]
#metabolome
second_view_pre=pre.iloc[:,784:918]
second_view_pre= pre_sample_id.join(second_view_pre)
#proteome
third_view_pre=pre.iloc[:,918:976]
third_view_pre=pre_sample_id.join(third_view_pre)
#gut
forth_view_pre=pre.iloc[:,976:983]
forth_view_pre=pre_sample_id.join(forth_view_pre)
#nares
fifth_view_pre=pre.iloc[:,983:]
fifth_view_pre=pre_sample_id.join(fifth_view_pre)
#---------------------------------------------------
def MofaData (x , y, view_name, group_name):
    ####df=my_data_withoutSampleID_into_Dataframeformat_dff=all_data_in_dataframe#####
    x=x.loc[:, x.columns != 'SampleID']
    #Here_write_the_hedear_name_of_sample_column_name
    column_count=len(x.columns)
    column_length=len(pd.DataFrame(x[[x.columns[1]]]))
    ########extract_value_column#########
    value_column=pd.concat([x, x.T.stack().reset_index(name='value')['value']], axis=1)
    value_column=value_column[['value']]
    value_column=pd.DataFrame(value_column)
    value_column_length=len(value_column[['value']])
    ########extract_view_and_group_column#########
    loop_value=0
    view_column=[]
    group_column=[]
    while loop_value < value_column_length:
        view_column.append(group_name)
        #here_you_can_insert_your_view_name
        group_column.append(view_name)
        #here_you_can_insert_your_group_name
        loop_value=loop_value+1
    view_column=pd.DataFrame(view_column, columns=["view"])
    group_column=pd.DataFrame(group_column, columns=["group"])
    ########extract_sample_column#########
    sample_column=[]
    for i in range(column_count):
        sample_column.append(y.loc[: ,'SampleID'])
    sample=pd.DataFrame(sample_column).T

    sample=pd.concat([sample, sample.T.stack().reset_index(name='sample')['sample']], axis=1)

    sample=sample[['sample']]
    sample_column=pd.DataFrame(sample)
    ########extract_feature_column#########
    feature_column=[]
    columns_name= x.columns
    for i in range(column_count):
        for _ in range(column_length):
            feature_column.append(columns_name[i+0]) 
    feature_column=pd.DataFrame(feature_column, columns=["feature"])
    ########concatonate_and_exporting_data#########
    concat=pd.concat([sample_column, group_column, feature_column, view_column, value_column],axis=1, sort=False)
    return concat 
#---------------------------------------------------
#Control_view
first_view_con = MofaData(first_view_con, first_view_con, 'Control', 'Deseq')
second_view_con = MofaData(second_view_con, second_view_con, 'Control', 'metabolome')
third_view_con = MofaData(third_view_con, third_view_con, 'Control', 'proteome')
forth_view_con = MofaData(forth_view_con, forth_view_con, 'Control', 'gut')
fifth_view_con = MofaData(fifth_view_con, fifth_view_con, 'Control', 'nares')
#Prediabetic_view
first_view_pre = MofaData(first_view_pre, first_view_pre, 'predipedic', 'Deseq')
second_view_pre = MofaData(second_view_pre, second_view_pre, 'predipedic', 'metabolome')
third_view_pre = MofaData(third_view_pre, third_view_pre, 'predipedic', 'proteome')
forth_view_pre = MofaData(forth_view_pre, forth_view_pre, 'predipedic', 'gut')
fifth_view_pre = MofaData(fifth_view_pre, fifth_view_pre, 'predipedic', 'nares')

conc=pd.concat([first_view_con, first_view_pre, second_view_con, second_view_pre, third_view_con, third_view_pre, forth_view_con, forth_view_pre, fifth_view_con, fifth_view_pre], axis=0)

