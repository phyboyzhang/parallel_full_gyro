ADD_LIBRARY(tp_preparation tp_preparation.F90)
Target_Link_Libraries(tp_preparation paradata_type picutilities piclayout constants parautilities 
                      outputdata)
