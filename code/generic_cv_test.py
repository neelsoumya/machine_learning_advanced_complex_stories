#####################################################
# Autoencoder using keras and tensorflow
#       applied to a dataset from the UCI machine
#       learning repository
#		NOTE: this takes in a mix of numeric and categorical features
#		numeric features should scaled before passing it to this program
#       NOTE: takes in many categoreical features and one numeric 
#			also has MLP on reduced features
#          this has class contrastive for deep learning
#       NOTE: this has more input parameters and is used to do model selection and print cross validation cv error cv loss
#
# INSTALLATION:
#   pip3 install tensorflow
#   https://www.tensorflow.org/install/install_mac
#   sudo pip3 install keras
#   https://keras.io/#installation
#   pip3 install graphviz
#   pip3 install pydot
#   pip3 install -U scikit-learn
#   pip3 install deap update_checker tqdm stopit
#   pip3 install tpot
#   pip3 install xgboost
#   pip3 install yellowbrick
#
#
# Usage:
#   python generic_cv_test.py df_med_suicide_features_lengthened_final_ml_f20_FULL_column_deleted.csv  drug_suidata_autoencoder_f20_MLP_FULL_column_deleted  0.5  relu  sigmoid
#  
# Adapted from:
#   https://blog.keras.io/building-autoencoders-in-keras.html
#   https://machinelearningmastery.com/display-deep-learning-model-training-history-in-keras/
#   https://machinelearningmastery.com/use-keras-deep-learning-models-scikit-learn-python/
#   https://github.com/rasbt/python-machine-learning-book/blob/master/code/ch06/ch06.ipynb
#
#####################################################


###################################################
# Load libraries
###################################################
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation
from keras.models import Model
from keras.utils import to_categorical
from keras import regularizers
from keras.optimizers import SGD
from keras.utils.vis_utils import plot_model
import keras
import pydot
import graphviz

import sklearn
# import scikit
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_val_score
from sklearn.metrics import roc_auc_score, accuracy_score, f1_score, matthews_corrcoef, average_precision_score
from sklearn.metrics import roc_curve, auc, precision_recall_fscore_support
from scipy import interp

import numpy as np
import pandas as pd
import pdb
import sys
import csv

import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split

from keras.layers import Input
from keras.datasets import mnist
from yellowbrick.text import TSNEVisualizer
from sklearn.ensemble import RandomForestClassifier
import tpot
from tpot import TPOTClassifier

from tkinter import *
from itertools import *

sns.set()

def function_tsne_generic(f_matrix, perplexities, str_save_file_tsne):

    """
    Function for tSNE
    Returns:

    Code adapted from:
    https://scikit-learn.org/stable/auto_examples/manifold/plot_t_sne_perplexity.html

    """

    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter
    from sklearn import manifold, datasets
    from yellowbrick.text import TSNEVisualizer

    # n_samples = 300
    n_components = 2
    i_num_panels_per_row = 5
    i_num_rows = int( np.ceil(len(perplexities)/i_num_panels_per_row) )

    (fig, subplots) = plt.subplots(i_num_rows, i_num_panels_per_row, figsize=(15, 8))

    # X, y = datasets.make_circles(n_samples=n_samples, factor=.5, noise=.05)

    # str_file_name = "breast-cancer-wisconsin_MOD.data"
    ## X_withcode = pd.read_csv(str_file_name)#, header=None)

    # modifications here and data munging
    # X = X_withcode.iloc[:,1:] # ignore code column
    X = f_matrix
    # pdb.set_trace()

    for i, perplexity in enumerate(perplexities):

        i_row_number = i%i_num_panels_per_row
        j_col_number = int( np.floor(i/i_num_panels_per_row) )

        #print(i)
        #print(j_col_number, i_row_number)

        ax = subplots[j_col_number][i_row_number]

        # t0 = time()
        tsne = manifold.TSNE(n_components=n_components, init='random',
                             random_state=0, perplexity=perplexity)
        Y = tsne.fit_transform(X)
        # t1 = time()
        # print("S-curve, perplexity=%d in %.2g sec" % (perplexity, t1 - t0))

        ax.set_title("Perplexity=%d" % perplexity)
        ax.scatter(Y[:, 0], Y[:, 1],  cmap=plt.cm.viridis) # c=color,
        ax.xaxis.set_major_formatter(NullFormatter())
        ax.yaxis.set_major_formatter(NullFormatter())
        ax.axis('tight')



    # plt.show()
    plt.savefig(str_save_file_tsne)
    plt.close()	

    # pdb.set_trace()

    return (Y)



def function_tsne_oneperplexity(f_matrix, perplexity, str_save_file_tsne):

    """
    Function for tSNE
    Returns:

    Code adapted from:
    https://scikit-learn.org/stable/auto_examples/manifold/plot_t_sne_perplexity.html

    """

    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter
    from sklearn import manifold, datasets
    from yellowbrick.text import TSNEVisualizer

    # n_samples = 300
    n_components = 2
    i_num_panels_per_row = 5
    i_num_rows = int( np.ceil(len(perplexities)/i_num_panels_per_row) )

    # (fig, subplots) = plt.subplots(i_num_rows, i_num_panels_per_row, figsize=(15, 8))

    # X, y = datasets.make_circles(n_samples=n_samples, factor=.5, noise=.05)

    # str_file_name = "breast-cancer-wisconsin_MOD.data"
    ## X_withcode = pd.read_csv(str_file_name)#, header=None)

    # modifications here and data munging
    # X = X_withcode.iloc[:,1:] # ignore code column
    X = f_matrix
    # pdb.set_trace()

    #for i, perplexity in enumerate(perplexities):

    #    i_row_number = i%i_num_panels_per_row
    #    j_col_number = int( np.floor(i/i_num_panels_per_row) )

    #    #print(i)
    #    #print(j_col_number, i_row_number)

    #    ax = subplots[j_col_number][i_row_number]

    #    # t0 = time()
    tsne = manifold.TSNE(n_components=n_components, init='random', random_state=0, perplexity=perplexity)
    Y = tsne.fit_transform(X)
    # t1 = time()
    # print("S-curve, perplexity=%d in %.2g sec" % (perplexity, t1 - t0))

    plt.figure()
    plt.grid()	
    plt.title("Perplexity=%d" % perplexity)
    plt.scatter(Y[:, 0], Y[:, 1],  cmap=plt.cm.viridis) # c=color,
    #plt.xaxis.set_major_formatter(NullFormatter())
    #plt.yaxis.set_major_formatter(NullFormatter())
    plt.axis('tight')



    # plt.show()
    plt.savefig(str_save_file_tsne)
    plt.close()	

    # pdb.set_trace()

    return (Y)



##############################################################
# functions
# 'btn_click' function continuously updates the input field
# whenever you enters a number
##############################################################
def btn_click(item):
    global expression
    expression = expression + str(item)
    print (item)

    # what is the new value in the tex field
    print ( input_text.get() )

    # set it to something if you want
    #input_text.set(expression)

    global i_temp_patient_index_1
    i_temp_patient_index_1 = input_text.get()

    # TODO: call AI function


# 'btn_clear' function clears the input field
def btn_clear():
    global expression
    expression = ""
    input_text.set("")



if __name__ == '__main__':

    ###################################################
    # Parse input command line arguments
    ###################################################
    str_input_filename         = str(sys.argv[1])   # 'df_drug_data_coded_tobesaved_curated.tsv'
    str_output_filename_prefix = str(sys.argv[2])   # 'drug_data_'
    f_dim_reduction_frac       = float(sys.argv[3]) # 0.4 # fraction by which dimensionality is to be reduced
    str_activation_layer1      = str(sys.argv[4])   # activation function for layer 1
    str_activation_layer2      = str(sys.argv[5])   # activation function for layer 2  
    
    ###################################################
    # Model parameters
    ###################################################
    # see https://machinelearningmastery.com/dropout-regularization-deep-learning-models-keras/
    f_split_train_test_fraction = 0.70 # 70% of data is for training
    # f_dim_reduction_frac = 0.4 # fraction by which dimensionality is to be reduced
    f_dropout                   = 0.20
    i_num_neurons_layer_1       = 10
    i_num_neurons_layer_2       = 4
    #str_activation_layer1       = 'relu'  # tanh
    #str_activation_layer2       = 'sigmoid'
    f_learning_rate             = 0.01
    f_learning_rate_decay       = 1e-6
    f_momentum                  = 0.9
    i_fitting_epochs            = 10000
    i_batch_size                = 320
    f_validation_split          = 0.33
    # k_kernel_initializer      = keras.initializers.VarianceScaling(scale=1.0, mode='fan_in', distribution='normal', seed=None)
    # k_kernel_initializer      = keras.initializers.glorot_uniform(seed=None) # Xavier Initialization
    k_kernel_initializer        = keras.initializers.Orthogonal(gain=1.0, seed=None)
    k_kernel_regularizer        = keras.regularizers.l1(0.00001)
    k_activity_regularizer      = keras.regularizers.l1(0.00001)
    # k_kernel_regularizer      = keras.regularizers.l1_l2(l1=0.01, l2=0.01)
    # sgd = keras.optimizers.SGD(lr=f_learning_rate, decay=f_learning_rate_decay,
    #                           momentum=f_momentum, nesterov=True)
    str_optimizer               = 'adadelta'
    str_loss_function           = 'kullback_leibler_divergence'#'mean_squared_error' #'binary_crossentropy' # 'kullback_leibler_divergence'
    # other loss functions
    # https://github.com/keras-team/keras/blob/master/keras/losses.py#L48


    #####################################################################################
    #  Load data
    #####################################################################################
    print("\n ********** Data Loading Section ********** \n")
    print("Loading dataset: \n")

    # read in pandas dataframe
    # generated from breast-cancer-wisconsin_MOD.data using data_munging.R
    temp_str_peptide_file = str_input_filename # "df_drug_data_coded_tobesaved_curated.tsv"  # "df_patients_joined_meds_orderedbymeds_curated_fewfields2.csv"
    temp_peptide_df = pd.read_csv(temp_str_peptide_file,  sep=",") # header=None,
    #            temp_peptide_df.columns = ['peptide', 'abundance']

    #####################################################################################
    # Split data into training and test set
    #####################################################################################
    print("\n ********** Train Test Split Section ********** \n")
    print("Splitting data into training and test set: \n")

    x_numeric       = temp_peptide_df.iloc[:, 0]  # temp_peptide_df["epithelial_cell_size"]
    x_categorical1  = temp_peptide_df.iloc[:, 1]
    x_categorical2  = temp_peptide_df.iloc[:, 2]  # temp_peptide_df["epithelial_cell_size"]
    x_categorical3  = temp_peptide_df.iloc[:, 3]
    x_categorical4  = temp_peptide_df.iloc[:, 4]
    x_categorical5  = temp_peptide_df.iloc[:, 5]
    x_categorical6  = temp_peptide_df.iloc[:, 6]
    x_categorical7  = temp_peptide_df.iloc[:, 7] 
    x_categorical8  = temp_peptide_df.iloc[:, 8]
    x_categorical9  = temp_peptide_df.iloc[:, 9]
    x_categorical10 = temp_peptide_df.iloc[:, 10]
    x_categorical11 = temp_peptide_df.iloc[:, 11]
    x_categorical12 = temp_peptide_df.iloc[:, 12]
    x_categorical13 = temp_peptide_df.iloc[:, 13]
    x_categorical14 = temp_peptide_df.iloc[:, 14]
    x_categorical15 = temp_peptide_df.iloc[:, 15]
    x_categorical16 = temp_peptide_df.iloc[:, 16]
    x_categorical17 = temp_peptide_df.iloc[:, 17]
    x_categorical18 = temp_peptide_df.iloc[:, 18]
    x_categorical19 = temp_peptide_df.iloc[:, 19]
    # x_categorical20 = temp_peptide_df.iloc[:, 20]
    y = temp_peptide_df.iloc[:, -1]  # temp_peptide_df["class"]

    x_categorical1 = keras.utils.to_categorical(x_categorical1)
    x_categorical2 = keras.utils.to_categorical(x_categorical2)
    x_categorical3 = keras.utils.to_categorical(x_categorical3) 
    x_categorical4 = keras.utils.to_categorical(x_categorical4) 
    x_categorical5 = keras.utils.to_categorical(x_categorical5)
    x_categorical6 = keras.utils.to_categorical(x_categorical6)
    x_categorical7 = keras.utils.to_categorical(x_categorical7)
    x_categorical8 = keras.utils.to_categorical(x_categorical8)
    x_categorical9 = keras.utils.to_categorical(x_categorical9)
    x_categorical10 = keras.utils.to_categorical(x_categorical10)
    x_categorical11 = keras.utils.to_categorical(x_categorical11)	
    x_categorical12 = keras.utils.to_categorical(x_categorical12)
    x_categorical13 = keras.utils.to_categorical(x_categorical13)
    x_categorical14 = keras.utils.to_categorical(x_categorical14)
    x_categorical15 = keras.utils.to_categorical(x_categorical15)
    x_categorical16 = keras.utils.to_categorical(x_categorical16)
    x_categorical17 = keras.utils.to_categorical(x_categorical17)
    x_categorical18 = keras.utils.to_categorical(x_categorical18)
    x_categorical19 = keras.utils.to_categorical(x_categorical19)
    # x_categorical20 = keras.utils.to_categorical(x_categorical20)	

    # x_concat = np.concatenate([x_numeric, x_categorical1, x_categorical2], axis = 1)

    ###############################################################################
    # TODO: concatenate all numeric and categorical features based on a metadata column
    # pdb.set_trace()
    x_concat = np.concatenate([np.transpose((np.vstack((x_numeric, x_numeric)))), x_categorical1, x_categorical2, x_categorical3, x_categorical4, x_categorical5, x_categorical6, x_categorical7, x_categorical8, x_categorical9, x_categorical10, x_categorical11, x_categorical12, x_categorical13, x_categorical14, x_categorical15, x_categorical16, x_categorical17, x_categorical18, x_categorical19], axis=1) # , x_categorical4

    # CAUTION: setting seed here 
    # TODO: remove seed later
    #import random	
    #random.seed(10)	
    x_train, x_test, y_train, y_test = train_test_split(x_concat,
                                                        y,
                                                        train_size = f_split_train_test_fraction,
                                                        test_size = 1 - f_split_train_test_fraction, random_state=10)

    # save these variables
    x_concat_orig = x_concat
    x_train_orig  = x_train
    x_test_orig   = x_test
    y_train_orig  = y_train
    y_test_orig   = y_test


    print(np.shape(x_test))
    print(np.shape(x_train))
    # pdb.set_trace()

    # find size of data and determine size of training set
    temp = np.array(temp_peptide_df)

    i_original_dim = np.shape(x_test)[1] # 700

    # save numpy file x_test for TESTING
    np.savetxt('x_test.csv', x_test, delimiter = ',')
	


    #####################################################################################
    # Feature scaling
    #####################################################################################
    print("\n ********** Data Munging Section ********** \n")
    print("Performing feature scaling: \n")

    #x_train_array = keras.utils.normalize(x_train_array)
    # feature scaling
    #x_test_array = keras.utils.normalize(x_test_array)
    # y_test = to_categorical(y_test)#,2)
    # x_train.shape
    # x_train.shape[1:] # 28 x 28
    # np.prod(x_train.shape[1:])


    #####################################################################################
    # Feature scaling
    #####################################################################################
    print("\n ********** Data Munging Section ********** \n")
    print("Performing feature scaling: \n")

    # flatten to vector for each image (60000 x 784)
    #tpl_flatten_new_dimensions = (  len(x_train),  np.prod(x_train.shape[1:])  )
    #x_train = np.reshape( x_train,  tpl_flatten_new_dimensions )

    #tpl_flatten_new_dimensions = (  len(x_test),  np.prod(x_test.shape[1:])  )
    #x_test = np.reshape( x_test,  tpl_flatten_new_dimensions )



    ###################################################
    # Initialize model
    ###################################################
    print("\n ********** Model Creation Section ********** \n")
    print("Creating autoencoder deep learning model: \n")

    i_encoding_dim = int( np.floor(f_dim_reduction_frac * i_original_dim) ) # 4

    input_img = Input(shape=(i_original_dim,))

    # goes from 8 to 4
    encoded = Dense(units=i_encoding_dim, activation=str_activation_layer1 ) #,
                    # kernel_regularizer=k_kernel_regularizer,
                    # activity_regularizer=k_activity_regularizer
                    # )

    # encoded object will take the input of input_img
    encoded(inputs=input_img)

    # back to number of dimensions you had originally
    # TODO: make parameyter
    # TODO: report cv error	
    decoded = Dense(units=i_original_dim, activation=str_activation_layer2)

    # decoded will take the output of encoder which is
    decoded(inputs=encoded(inputs=input_img))

    ######################################################################################################
    # these are all inputs

    # we need to create the model
    #   first is input of class Input
    #   second is output Tensor (not of class Input)
    #
    #   CONCEPT: when you create the model, you define the path that the data will follow,
    #       which in this case is from the input to the output:
    #   See:
    #       https://stackoverflow.com/questions/43765381/anyone-some-api-about-the-keras-layers-input
    ######################################################################################################
    model_autoencoder = Model( input_img, decoded(inputs=encoded(inputs=input_img)) )

    model_autoencoder.compile(loss=str_loss_function,
                              optimizer=str_optimizer,
                              metrics=['accuracy'])

    model_autoencoder.summary()

    ###################################################
    # Now create decoder model
    # Before you create that recall that the input to decoder model
    #       needs to be of class Input

    # so create an Input class for decoder
    ###################################################

    input_to_decoder = Input(shape=(i_encoding_dim,))


    model_autoencoder.layers

    # Last layer of autoencoder Dense object
    dense_decoder = model_autoencoder.layers[-1]

    # NOTE: dense_decoder is an object of type Dense()
    # hence you need to give it an object of type Input()
    # In this case this will be the output of encoder
    #       also called input_to_decoder

    dense_decoder(input_to_decoder)

    # Putting it altogether

    # Now decoder model has
    #   input:  input_to_decoder
    #   output:
    model_decoder = Model(input_to_decoder, dense_decoder(input_to_decoder)  )



    #####################################################################################
    # Plot model fit over training epochs
    #  https://machinelearningmastery.com/display-deep-learning-model-training-history-in-keras/
    #####################################################################################

    # Fit model
    # TODO: print cv error
    #    https://machinelearningmastery.com/evaluate-performance-deep-learning-models-keras/	
    #    https://datascience.stackexchange.com/questions/25267/keras-difference-beetween-val-loss-and-loss-during-training
    history = model_autoencoder.fit(x_train,
                                    x_train,
                                    validation_split=f_validation_split,
                                    epochs=i_fitting_epochs,
                                    batch_size=i_batch_size,
                                    verbose=0)
    print('Cross validation error ...... \n')
    print( history.history['val_loss'][i_fitting_epochs-1] )	
    print(' ....   \n')
    np.savetxt(str_output_filename_prefix + '_xval_loss_test.csv', [ history.history['val_loss'][i_fitting_epochs-1] ], delimiter = ',')	

    # list all data in history
    # print(history.history.keys())
    # summarize history for accuracy
    # plt.figure()
    # plt.grid()
    # plt.plot(history.history['acc'])
    # plt.plot(history.history['val_acc'])
    # plt.title('model accuracy')
    # plt.ylabel('accuracy')
    # plt.xlabel('epoch')
    # plt.legend(['train', 'test'], loc='upper left')
    # plt.tight_layout()
    # plt.savefig(str_output_filename_prefix + 'autoencoder_learning_curve_accuracy.png', dpi=300)
    # plt.show()
    # summarize history for loss
    # plt.figure()
    # plt.grid()
    # plt.plot(history.history['loss'])
    # plt.plot(history.history['val_loss'])
    # plt.title('model loss')
    # plt.ylabel('loss')
    # plt.xlabel('epoch')
    # plt.legend(['train', 'test'], loc='upper left')
    # plt.tight_layout()
    # plt.savefig(str_output_filename_prefix + 'autoencoder_learning_curve_loss.png', dpi=300)
    # plt.show()


    #####################################################################################
    # Visualize encoder and decoder model output
    #       predict on test sets
    #   With a difference for autoencoders
    #   since we are not interested in the final model output
    #       just in the intermediate model output
    #####################################################################################

    #####################################################################################
    # This is also Generate predictions on new data:
    #####################################################################################
    print("\n ********** Model Prediction and Visualization Section ********** \n")
    print("Printing model prediction details on test set: \n")
    print("Visualize encoder and decoder model output on test set: \n")

    # predict on test set using encoder
    #  NOT autoencoder
    # the following is wrong
    # predictions_test_encoded = model_autoencoder.predict(x_test)

    # CONCEPT: the autoencoder model is the complete model consisting of
    #           BOTH encoder and decoder

    # we will now build an encoder Model
    #   input:  input_img
    #   output: encoded Dense object
    model_encoder = Model(input_img, encoded(inputs=input_img))

    # model_encoder.compile()
    predictions_testset_model_encoder = model_encoder.predict(x_test)

    # plot and visualize encoded matrix in tSNE space
    perplexities = [1, 2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120]
    tsne_reduced_matrix = function_tsne_generic(f_matrix=predictions_testset_model_encoder,
                          perplexities=perplexities,
                          str_save_file_tsne=str_output_filename_prefix + 'autoencoder_tsne_visualization.eps')
    # 'autoencoder_tsne_visualization.eps'

    # TESTING code plot and visualize encoded matrix in tSNE space
    perplexities = [20, 20, 20, 20, 20, 20, 20, 20]
    tsne_reduced_matrix_TESTING = function_tsne_generic(f_matrix=predictions_testset_model_encoder,
                          perplexities=perplexities,
                          str_save_file_tsne=str_output_filename_prefix + 'autoencoder_tsne_visualization_TESTING.png')


    tsne_with_labels = TSNEVisualizer()
    #tsne_with_labels.fit(predictions_testset_model_encoder,y_test)
    #tsne_with_labels.poof()

    # save this in a file also
    f_file_save = open(str_output_filename_prefix + '_tsne_save.csv', 'w')
    for item in predictions_testset_model_encoder:
        f_file_save.writelines(str(item))
        f_file_save.writelines('\n')

    f_file_save.close()

    f_file_save = open(str_output_filename_prefix + '_tsne_reduced_matrix_save.csv', 'w')
    for item in tsne_reduced_matrix:
        f_file_save.writelines(str(item))
        f_file_save.writelines('\n')

    f_file_save.close()

    # plot scatter plot of these tsne latent dimensions or encoder weights on test set labels		
    plt.figure(figsize=(7,5))	
    plt.scatter(tsne_reduced_matrix[:,0], tsne_reduced_matrix[:,1], c=y_test)
    plt.colorbar()
    plt.savefig(str_output_filename_prefix + '_tsne_reduced_matrix_colored_by_label.png', dpi=300)
	
    # TODO: Visualize and predict tsne clusters based on features
    #  use x_test and tsne_reduced_matrix_TESTING USE np.concatenate	axis = 1

    # conversational AI
    b_conversation = True
    i_temp_perplexity      = 20
    i_temp_patient_index_1 = 34
    # TODO: have an index patient against which tSNE co-ordinate change can be compared	
    i_temp_patient_index_INDEX = 100	
    i_temp_col_mutate_1    = 6
    i_temp_col_mutate_2    = 8


    # TODO: call AI function
    x_test_archive = x_test
	
    # TODO: check if 1 before setting to 0 i.e complement	


    # TODO: repeat this for MLP simple deep learning prediction
    while b_conversation:
        # pdb.set_trace()
        # adapt or mutate
        x_test[i_temp_patient_index_1,i_temp_col_mutate_1] = 1 # mutate this patient to a 1
        x_test[i_temp_patient_index_1,i_temp_col_mutate_2] = 1 # mutate this patient to a 1
		
        # predict again using encoder
        predictions_testset_model_encoder = model_encoder.predict(x_test)
        # run tsne again, bear in mind this is stochastic and perturbing one point may change others
        perplexities = [i_temp_perplexity, i_temp_perplexity, i_temp_perplexity, i_temp_perplexity, i_temp_perplexity, i_temp_perplexity, i_temp_perplexity, i_temp_perplexity]
        tsne_reduced_matrix_TESTING = function_tsne_generic(f_matrix=predictions_testset_model_encoder, perplexities=perplexities, str_save_file_tsne=str_output_filename_prefix + 'autoencoder_tsne_visualization_TESTING2.png')
        tsne_reduced_matrix_TESTING2 = function_tsne_oneperplexity(f_matrix=predictions_testset_model_encoder, perplexity=i_temp_perplexity, str_save_file_tsne=str_output_filename_prefix + 'autoencoder_tsne_visualization_TESTING2.png')
	
        # does this change the prediction for 34th patient
        tsne_reduced_matrix_TESTING[i_temp_patient_index_1,]
        tsne_reduced_matrix_TESTING2[i_temp_patient_index_1,]
        print(tsne_reduced_matrix_TESTING2[i_temp_patient_index_1,])
        # TODO: logic for co-ordinates remember that distancee does not matter in tsne space
        # but still have logic if points are "distinct" suggest looking at which quadrant the points are in
        # and also include printing index patient	 	
	
        # printing INDEX patient 		
        print(tsne_reduced_matrix_TESTING2[i_temp_patient_index_INDEX,])		

        # perturb back
        x_test[i_temp_patient_index_1,i_temp_col_mutate_1] = 0 # mutate this patient to a 1
        x_test[i_temp_patient_index_1,i_temp_col_mutate_2] = 0 # mutate this patient to a 1
		
        # predict again using encoder
        predictions_testset_model_encoder = model_encoder.predict(x_test)
        # run tsne again, bear in mind this is stochastic and perturbing one point may change others
        tsne_reduced_matrix_TESTING2 = function_tsne_oneperplexity(f_matrix=predictions_testset_model_encoder, perplexity=i_temp_perplexity, str_save_file_tsne=str_output_filename_prefix + 'autoencoder_tsne_visualization_TESTING2_back.png')
        print(tsne_reduced_matrix_TESTING2[i_temp_patient_index_1,])
		
        # printing INDEX patient 		
        print(tsne_reduced_matrix_TESTING2[i_temp_patient_index_INDEX,])



        #####################
        # create a window
        #####################
        #window = Tk()
        #window.geometry("312x324")
        #window.title("Conversational AI")

        # text box and frame
        #expression = ""
        # 'StringVar()' is used to get the instance of input field
        #input_text = StringVar()

        # creating a frame for the input field
        #input_frame = Frame(window, width=312, height=50, bd=0, highlightbackground="black",
        #                    highlightcolor="black", highlightthickness=1)
        #input_frame.pack(side=TOP)

        # creating a input field inside the 'Frame'
        #input_field = Entry(input_frame, font=('arial', 18, 'bold'), textvariable=input_text,
        #                    width=50, bg="#eee", bd=0, justify=RIGHT)
        #input_field.grid(row=0, column=0)
        #input_field.pack(ipady=10)  # 'ipady' is internal padding to increase the height of input field

        # creating another 'Frame' for the button below the 'input_frame'
        #btns_frame = Frame(window, width=312, height=272.5, bg="grey")
        #btns_frame.pack()

        # first row
        #divide = Button(btns_frame, text="Compute", fg="black", width=10, height=3, bd=0, bg="#eee",
        #                cursor="hand2", command=lambda: btn_click("123")).grid(row=0, column=3, padx=1, pady=1)

        # clear = Button(btns_frame, text="Clear", fg="black", width=32, height=3, bd=0,
        #               bg="#eee", cursor="hand2",
        #               command=lambda: btn_clear()).grid(row=0, column=0, columnspan=3, padx=1, pady=1)

        #window.mainloop()

        #print (input_field)

        # TODO: restore x_test
        # x_test = x_test_archive		

        pdb.set_trace()



    # TODO: natural language conversation if age different then not die
    # TODO: conversation with itself and explor the whole feature space and create a story
    #           broad overarching narrative and stopry that if you change this featie and that feature
    #           then for some patients this happens
	
    # TODO: in R then use class contrastive OR edit distance based on how many features changed	
    # TODO: try combined loss functions KL + SSR
    # TODO: use edit distance based metric for graph deep learning
    # TODO:      if edit dustnace is based on tree distance, then maybe better than plain autoencoder
    #            sicne say 4 edits is wha t defjes clusters in tsne then on tree ditance will be less

    # TESTING code
    f_file_save = open(str_output_filename_prefix + '_tsne_reduced_matrix_save_TESTING.csv', 'w')
    for item in tsne_reduced_matrix_TESTING:
        f_file_save.writelines(str(item))
        f_file_save.writelines('\n')

    f_file_save.close()

    # plot scatter plot of these tsne latent dimensions or encoder weights on test set labels		
    # TODO: also color by each feature say an important feature like dementia	
    plt.figure(figsize=(7,5))	
    plt.scatter(tsne_reduced_matrix_TESTING[:,0], tsne_reduced_matrix_TESTING[:,1], c=y_test)
    plt.colorbar()
    plt.savefig(str_output_filename_prefix + '_tsne_reduced_matrix_colored_by_label_TESTING.png', dpi=300)
	
	
    # plot scatter plot of these latent dimensions or encoder weights on test set
    plt.figure(figsize=(7, 5)) 
    plt.scatter(predictions_testset_model_encoder[:,0], predictions_testset_model_encoder[:,1], c=y_test)
    plt.colorbar()
    plt.savefig(str_output_filename_prefix + '_reduced_dimensions_plot.png', dpi=300)
 		

    # feed these predictions to decoder
    predictions_testset_model_decoder = model_decoder.predict(predictions_testset_model_encoder)




    #####################################################################################
    # Visualize balance or imbalance of training data
    #####################################################################################
    plt.figure(figsize=(8, 4))
    sns.countplot(x=y_train)
    plt.savefig(str_output_filename_prefix + '_autoencoder_balance_trainingset.png', dpi=300)

    plt.figure(figsize=(8, 4))
    sns.countplot(x=y_test)
    plt.savefig(str_output_filename_prefix + '_autoencoder_balance_testset.png', dpi=300)

    #####################################################################################
    # Print summary of a model and inspect model
    #####################################################################################
    print("\n ********** Model Summary Section ********** \n")
    print("Printing autoencoder model summary and model details: \n")
    print(model_autoencoder.summary())
    keras.utils.print_summary(model_autoencoder, line_length=None, positions=None, print_fn=None)
    print(model_autoencoder.output_shape)
    print(model_autoencoder.input_shape)
    print(model_autoencoder.get_config())
    print(model_autoencoder.get_weights())

    print("\n ********** Model Parameters ********** \n")
    print("Training and test split", f_split_train_test_fraction)
    print("Learning rate: ", f_learning_rate)
    print("Learning rate decay: ", f_learning_rate_decay)
    print("Optimizer:", "SGD")
    print("Momentum: ", f_momentum)
    print("Fitting epochs: ", i_fitting_epochs)
    print("Batch size: ", i_batch_size)
    print("Validation split of training set: ", f_validation_split)
    print("Dropout probability: ", f_dropout)
    print("Number of neurons in first hidden layer: ", i_num_neurons_layer_1)
    print("Number of neurons in second hidden layer: ", i_num_neurons_layer_2)
    print("Activation function for first hidden layer: ", str_activation_layer1)
    print("Activation function for first hidden layer: ", str_activation_layer2)
    print(
    "Kernel initialization: ", "Orthogonal")  # k_kernel_initializer   = keras.initializers.Orthogonal(gain=1.0, seed=None)
    print("Kernel regularization: ", "L1")  # k_kernel_regularizer   = keras.regularizers.l1(0.01)
    print("Kernel activity initialization: ", "L1")  # k_activity_regularizer = keras.regularizers.l1(0.01)
    print("Loss function: ", str_loss_function)  # k_activity_regularizer = keras.regularizers.l1(0.01)	
    print("**************************************\n")

    #####################################################################################
    # Plot model
    #####################################################################################
    # plot_model(model, to_file='model.png')
    #keras.utils.plot_model(model_autoencoder,
    #                       to_file=str_output_filename_prefix + 'autoencoder_keras_model.png',
    #                       show_shapes=True,
    #                       show_layer_names=True,
    #                       rankdir='TB')




    #####################################################################################
    #
    #####################################################################################
	
	
    #####################################################################################
    # MLP on reduced features
    #####################################################################################

    ###################################################
    # Model parameters
    ###################################################
    # see https://machinelearningmastery.com/dropout-regularization-deep-learning-models-keras/
    # f_split_train_test_fraction = 0.75 # 75% of data is for training
    f_dropout                   = 0.50
    i_num_neurons_layer_1       = 4
    i_num_neurons_layer_2       = 2
    # str_activation_layer1       = 'relu'  # tanh
    # str_activation_layer2       = 'sigmoid'
    f_learning_rate             = 0.01
    f_learning_rate_decay       = 1e-6
    f_momentum                  = 0.9
    i_fitting_epochs            = 10000
    i_batch_size                = 500
    f_validation_split          = 0.33
    # k_kernel_initializer      = keras.initializers.VarianceScaling(scale=1.0, mode='fan_in', distribution='normal', seed=None)
    # k_kernel_initializer      = keras.initializers.glorot_uniform(seed=None) # Xavier Initialization
    k_kernel_initializer        = keras.initializers.Orthogonal(gain=1.0, seed=None)
    k_kernel_regularizer        = keras.regularizers.l1(0.01)
    k_activity_regularizer      = keras.regularizers.l1(0.01)
    # k_kernel_regularizer      = keras.regularizers.l1_l2(l1=0.01, l2=0.01)
    sgd = keras.optimizers.SGD(lr=f_learning_rate, decay=f_learning_rate_decay,
                           momentum=f_momentum, nesterov=True)
	
	
    # feed these predictions of encoder model on all data
    predictions_ALL_model_encoder = model_encoder.predict(x_concat)
    x_concat = predictions_ALL_model_encoder
    x_train, x_test, y_train, y_test = train_test_split(x_concat,
                                                    y,
                                                    train_size=f_split_train_test_fraction,
                                                    test_size=1 - f_split_train_test_fraction)

    # x_train = model_encoder.predict(x_train)
    # x_test  = model_encoder.predict(x_test)
    # y_train = y_train
    # y_test  = y_test


    print(np.shape(x_train))
    print(np.shape(x_test))

    #####################################################################################
    # Feature scaling
    #####################################################################################
    print("\n ********** Data Munging Section ********** \n")
    print("Performing feature scaling: \n")


    ###################################################
    # Initialize model
    ###################################################
    print("\n ********** Model Creation Section ********** \n")
    print("Creating multi layer perceptron deep learning model: \n")

    model = Sequential()

    ###################################################
    # Adding layers
    # Stacking layers is as easy as .add():
    ###################################################
    # model.add(Dense(units=64, activation='relu', input_dim=9))
    i_input_dimension = np.shape(x_test)[1]
    model.add(Dense(units=i_num_neurons_layer_1, activation=str_activation_layer2, input_dim=i_input_dimension,
                kernel_initializer=k_kernel_initializer,
                kernel_regularizer=k_kernel_regularizer,
                activity_regularizer=k_activity_regularizer
                )
          )
    #model.add(Dropout(f_dropout))
    #model.add(Dense(units=i_num_neurons_layer_2, activation=str_activation_layer2))


    model.compile(loss=keras.losses.sparse_categorical_crossentropy,
              optimizer=sgd,
              metrics=['accuracy'])

    # model.compile(loss=keras.losses.categorical_crossentropy,
    #              optimizer=keras.optimizers.SGD(lr=0.01, momentum=0.9, nesterov=True))



    #####################################################################################
    # Plot model fit over training epochs
    #  https://machinelearningmastery.com/display-deep-learning-model-training-history-in-keras/
    #####################################################################################

    # Account for imbalanced classes in model fitting
    # NOTE: Currently only works for binary classification with labels 1 and 0 (numeric)

    i_class_weight_label_1 = sum([int(x) for x in y_train]) / len(y_train) # number of 1s in target (y) of training set
    i_class_weight_label_0 = 1 - i_class_weight_label_1
    dict_class_weights     = {1:i_class_weight_label_1, 0:i_class_weight_label_0} # dict of class weights

    history = model.fit(x=x_train, y=y_train,
                    validation_split=f_validation_split,
                    class_weight=dict_class_weights,
                    epochs=i_fitting_epochs,
                    batch_size=i_batch_size,
                    verbose=0)

    # list all data in history
    # print(history.history.keys())
    # summarize history for accuracy
    # plt.figure()
    # plt.grid()
    # plt.plot(history.history['acc'])
    # plt.plot(history.history['val_acc'])
    # plt.title('model accuracy')
    # plt.ylabel('accuracy')
    # plt.xlabel('epoch')
    # plt.legend(['train', 'test'], loc='upper left')
    # plt.tight_layout()
    # plt.savefig(str_output_filename_prefix + 'multipledatatypes_learning_curve_accuracy_drugmedsuicide_MLP_f20FULL.png', dpi=300)
    # plt.show()

    # summarize history for loss
    # plt.figure()
    # plt.grid()
    # plt.plot(history.history['loss'])
    # plt.plot(history.history['val_loss'])
    # plt.title('model loss')
    # plt.ylabel('loss')
    # plt.xlabel('epoch')
    # plt.legend(['train', 'test'], loc='upper left')
    # plt.tight_layout()
    # plt.savefig(str_output_filename_prefix + 'multipledatatypes_learning_curve_loss_drugmedsuicide_MLP_f20FULL.png', dpi=300)
    # plt.show()


    #####################################################################################
    # Evaluate your performance in one line:
    #####################################################################################

    print("\n ********** Model Evaluation Section ********** \n")
    print("Printing model performance details on test set: \n")

    loss_and_metrics = model.evaluate(x_test, y_test, batch_size=i_batch_size)
    print("loss_and_metrics", loss_and_metrics)

    #####################################################################################
    # TODO: MODEL SELECTION
    # https://machinelearningmastery.com/use-keras-deep-learning-models-scikit-learn-python/
    #####################################################################################

    #####################################################################################
    # Generate predictions on new data:
    #####################################################################################
    print("\n ********** Model Prediction Section ********** \n")
    print("Printing model prediction details on test set: \n")

    classes = model.predict(x_test, batch_size=i_batch_size)
    print(y_test)
    print(classes)
    model.predict_classes(x_test, batch_size=i_batch_size, verbose=0)  # , steps=1)
    # pdb.set_trace()

    #####################################################################################
    # Plot AUPR curves or ROC curves
    # https://scikit-learn.org/stable/auto_examples/model_selection/plot_precision_recall.html#sphx-glr-auto-examples-model-selection-plot-precision-recall-py	
    # https://github.com/rasbt/python-machine-learning-book/blob/master/code/ch06/ch06.ipynb
    #####################################################################################
    # TODO: AUPR curves not ROC
    # see https://hackernoon.com/simple-guide-on-how-to-generate-roc-plot-for-keras-classifier-2ecc6c73115a
    probas = model.predict(x_test)  # x_test

    # pdb.set_trace()
    # fpr, tpr, thresholds = roc_curve(y_test[:,1], probas[:,1], pos_label=1)
    fpr, tpr, thresholds = roc_curve(y_test, probas[:, 1], pos_label=1)

    mean_tpr = 0.0
    mean_fpr = np.linspace(0, 1, 100)
    all_tpr = []

    mean_tpr += interp(mean_fpr, fpr, tpr)
    mean_tpr[0] = 0.0
    roc_auc = auc(fpr, tpr)

    print("ROC score", roc_auc)

    plt.figure(figsize=(7, 5))
    plt.grid()
    plt.plot(fpr, tpr, lw=1, label='ROC (area = %0.2f)' % (roc_auc))
    plt.plot([0, 1], [0, 1], linestyle='--', color=(0.6, 0.6, 0.6), label='random guessing')
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    plt.title('Receiver Operator Characteristic')
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(str_output_filename_prefix + 'multipledatatypes_roc_curve_drugmedsuicide_f20FULL_MLP.png', dpi=300)



    #####################################################################################
    # Visualize balance or imbalance of training data
    #####################################################################################
    plt.figure(figsize=(8, 4))
    sns.countplot(x=y_train)
    plt.savefig(str_output_filename_prefix + 'multipledatatypes_balance_trainingset_drugmedsuicide_f20FULL_MLP.png', dpi=300)

    plt.figure(figsize=(8, 4))
    sns.countplot(x=y_test)
    plt.savefig(str_output_filename_prefix + 'multipledatatypes_balance_testset_drugmedsuicide_f20FULL_MLP.png', dpi=300)

    #####################################################################################
    # Print summary of a model and inspect model
    #####################################################################################
    print("\n ********** Model Summary Section ********** \n")
    print("Printing feedforward neural network built on autoencoder model summary and model details: \n")
    print(model.summary())
    keras.utils.print_summary(model, line_length=None, positions=None, print_fn=None)
    print(model.output_shape)
    print(model.input_shape)
    print(model.get_config())
    #print(model.get_weights())

    print("\n ********** Model Parameters ********** \n")
    print("\nTraining and test split", f_split_train_test_fraction)
    print("Learning rate: ", f_learning_rate)
    print("Learning rate decay: ", f_learning_rate_decay)
    print("Optimizer:", "SGD")
    print("Momentum: ", f_momentum)
    print("Fitting epochs: ", i_fitting_epochs)
    print("Batch size: ", i_batch_size)
    print("Validation split of training set: ", f_validation_split)
    print("Dropout probability: ", f_dropout)
    print("Number of neurons in first hidden layer: ", i_num_neurons_layer_1)
    print("Number of neurons in second hidden layer: ", i_num_neurons_layer_2)
    print("Activation function for first hidden layer: ", str_activation_layer1)
    print("Activation function for first hidden layer: ", str_activation_layer2)
    print("Kernel initialization: ", "Orthogonal")  # k_kernel_initializer   = keras.initializers.Orthogonal(gain=1.0, seed=None)
    print("Kernel regularization: ", "L1")  # k_kernel_regularizer   = keras.regularizers.l1(0.01)
    print("Kernel activity initialization: ", "L1")  # k_activity_regularizer = keras.regularizers.l1(0.01)
    print("Loss function:" , "sparse_categorical_cross_entropy")
    print("**************************************")

    #####################################################################################
    # Plot model
    #####################################################################################
    # plot_model(model, to_file='model.png')
    #keras.utils.plot_model(model, to_file='multiple_datatypes_keras_model.png', show_shapes=True,
    #                       show_layer_names=True, rankdir='TB')




    #####################################################################################
    # Save model
    #####################################################################################
    print("\n ********** Model Save Section ********** \n")
    print("       Saving model ......  \n")
    str_model_save_filename = str_output_filename_prefix + "model_file_uci_dataset_drugmedsuicide_f20FULL_MLP.h5"
    model.save(str_model_save_filename)
    # model_saved = load_model(str_model_save_filename)

    ##################################################################
    # Other baseline models
    ##################################################################

    # Random forest
    rf = RandomForestClassifier(max_depth=3, n_estimators=10)
    rf.fit(x_train, y_train)

    y_pred_rf = rf.predict_proba(x_test)[:, 1]
    fpr_rf, tpr_rf, thresholds_rf = roc_curve(y_test, y_pred_rf)
    auc_rf = auc(fpr_rf, tpr_rf)
    print("AUC from a baseline random forest algorithm: ", auc_rf)
    pdb.set_trace()
    # f1_score_rf = f1_score(y_test, y_pred_rf)
    # print("F1 score from a baseline random forest algorithm: ", f1_score_rf)
    # matthews_corrcoef_rf = matthews_corrcoef(y_test, y_pred_rf)
    # print("Matthews correlation coefficient from a baseline random forest algorithm: ", matthews_corrcoef_rf)
    # precision, recall, fbeta_score, support = precision_recall_fscore_support(y_test, y_pred_rf) #, average='weighted')
    # https://en.wikipedia.org/wiki/F1_score	
    # https://scikit-learn.org/stable/modules/generated/sklearn.metrics.average_precision_score.html
    # https://stats.stackexchange.com/questions/157012/area-under-precision-recall-curve-auc-of-pr-curve-and-average-precision-ap
    # https://scikit-learn.org/stable/modules/generated/sklearn.metrics.matthews_corrcoef.html
    # https://scikit-learn.org/stable/modules/generated/sklearn.metrics.f1_score.html
    # print("Precision recall F1 score from a baseline random forest algorithm: ", fbeta_score) 	


    ###################################################################
    # Class-contrastive for machine learning models
    ###################################################################     
    #  start for loop
    # TODO: 1. can mutate features in x_concat do for each feature
    # TODO: call AI function
    
    # save or archive original arrtay before mutating it
    x_test_original_archive = x_test_orig
	
    # 2 onwards are categorical columns
    # this is the offset, al categorical columns must be on or after this
    i_categorical_offset = 2
    (i_num_patients_testset, i_num_columns_categorical_testset) = np.shape(x_test_orig[:,i_categorical_offset:])
    
    # TODO: check if 1 before setting to 0 i.e complement	
    # adapt or mutate
    # i_temp_patient_index_1 = 10
    # i_temp_col_mutate_1    = 3
    
    # the number of features counted so far are one-hot encoded so divide by 2 since each feature is only 0/1
    #  and pytjon will conevrt it to 10 or 01
    i_num_columns_categorical_testset_absolute = int(i_num_columns_categorical_testset/2)
    
    # this is a data frame that will store all class-contrastive probabilities for each patient (row) and each categorical feature (column)
    df_class_contrastive_ml = np.zeros((i_num_patients_testset, i_num_columns_categorical_testset_absolute ))
    
    # for each patient
    for i_temp_patient_index_1 in range(0, i_num_patients_testset - 1):
        # for each categorical column
        for i_temp_col_mutate_1 in range(i_categorical_offset, round((i_num_columns_categorical_testset + i_categorical_offset - 1)/2) ):      
            # NOTE: do this for the original x_test that is x_test_orig
            #  each feature is one hot encoded adn each feature is binary so python will convert to 10 or 01
            x_test_orig[i_temp_patient_index_1,i_temp_col_mutate_1]   = 1 # mutate this patient to a 1
            x_test_orig[i_temp_patient_index_1,i_temp_col_mutate_1+1] = 0 # mutate this patient to a 1
    
            # reduce dimensions using autoencoder
            predictions_ALL_model_encoder_test = model_encoder.predict(x_test_orig)
    
            #  use a random forest on reduced dimensions
            y_pred_rf_contrast = rf.predict_proba(predictions_ALL_model_encoder_test)[:, 1]
            # TODO: 2. pick one patient and this is predicted probability
            y_pred_rf_contrast[i_temp_patient_index_1]
    
    
            # mutate back
            #  each feature is one hot encoded adn each feature is binary so python will convert to 10 or 01
            x_test_orig[i_temp_patient_index_1,i_temp_col_mutate_1]   = 0 # mutate this patient to a 1
            x_test_orig[i_temp_patient_index_1,i_temp_col_mutate_1+1] = 1 # mutate this patient to a 1
    
            # reduce dimensions using autoencoder
            predictions_ALL_model_encoder_test = model_encoder.predict(x_test_orig)
    
            #  use a random forest on reduced dimensions
            y_pred_rf_contrast_back = rf.predict_proba(predictions_ALL_model_encoder_test)[:, 1]
            # TODO: 2. pick one patient and this is predicted probability
            y_pred_rf_contrast_back[i_temp_patient_index_1]
    
            i_temp_delta_prob = y_pred_rf_contrast[i_temp_patient_index_1] - y_pred_rf_contrast_back[i_temp_patient_index_1]
    
            # store i_temp_delta_prob in a data frame   
            #       data frame row is patient and column is feature
            # pdb.set_trace()
            df_class_contrastive_ml[i_temp_patient_index_1, i_temp_col_mutate_1 - i_categorical_offset] = i_temp_delta_prob        

    #  end for loop
    
    # restore test set
    x_test_orig = x_test_original_archive

    pdb.set_trace()

    # save data frame of class-contrastive results to hard disk
    # this will be read from R
    np.savetxt('df_class_contrastive_ml_from_python.csv', df_class_contrastive_ml, delimiter='\t')

    # list of all array indices (categorical only)
    list_temp = [x for x in range(i_categorical_offset, round((i_num_columns_categorical_testset + i_categorical_offset - 1) / 2)) ]
    # tuple of all combinations (2 at a time)
    tuple_temp = [x for x in combinations(list_temp, 2) ]

    with open("list_tuples.csv", 'w', newline='\n') as outfile:
        wr = csv.writer(outfile)#, quoting=csv.QUOTE_ALL)
        for (x, y) in tuple_temp:
            print(x,y)
            wr.writerow([x,y])
            
    outfile.close()        

    # this is a data frame that will store all class-contrastive probabilities for each patient (row) and each categorical feature (column) ALL combinations
    i_num_combinations_columns = len(tuple_temp)
    df_class_contrastive_ml_2comb = np.zeros((i_num_patients_testset, i_num_combinations_columns))

    # for each patient
    for i_temp_patient_index_1 in range(0, i_num_patients_testset - 1):
        # for each combination of categorical column
        
        #  set up a temporary row counter
        i_temp_row_counter = 0        
        
        for (i_temp_col_mutate_1, i_temp_col_mutate_2) in tuple_temp:
            # NOTE: do this for the original x_test that is x_test_orig
            #  each feature is one hot encoded adn each feature is binary so python will convert to 10 or 01
            x_test_orig[i_temp_patient_index_1, i_temp_col_mutate_1] = 1  # mutate this patient to a 1
            x_test_orig[i_temp_patient_index_1, i_temp_col_mutate_1 + 1] = 0  # mutate this patient to a 1

            x_test_orig[i_temp_patient_index_1, i_temp_col_mutate_2] = 1  # mutate this patient to a 1
            x_test_orig[i_temp_patient_index_1, i_temp_col_mutate_2 + 1] = 0  # mutate this patient to a 1

            # reduce dimensions using autoencoder
            predictions_ALL_model_encoder_test = model_encoder.predict(x_test_orig)

            #  use a random forest on reduced dimensions
            y_pred_rf_contrast = rf.predict_proba(predictions_ALL_model_encoder_test)[:, 1]
            # TODO: 2. pick one patient and this is predicted probability
            y_pred_rf_contrast[i_temp_patient_index_1]

            # mutate back
            #  each feature is one hot encoded adn each feature is binary so python will convert to 10 or 01
            x_test_orig[i_temp_patient_index_1, i_temp_col_mutate_1] = 0  # mutate this patient to a 1
            x_test_orig[i_temp_patient_index_1, i_temp_col_mutate_1 + 1] = 1  # mutate this patient to a 1

            x_test_orig[i_temp_patient_index_1, i_temp_col_mutate_2] = 0  # mutate this patient to a 1
            x_test_orig[i_temp_patient_index_1, i_temp_col_mutate_2 + 1] = 1  # mutate this patient to a 1

            # reduce dimensions using autoencoder
            predictions_ALL_model_encoder_test = model_encoder.predict(x_test_orig)

            #  use a random forest on reduced dimensions
            y_pred_rf_contrast_back = rf.predict_proba(predictions_ALL_model_encoder_test)[:, 1]
            # TODO: 2. pick one patient and this is predicted probability
            y_pred_rf_contrast_back[i_temp_patient_index_1]

            i_temp_delta_prob = y_pred_rf_contrast[i_temp_patient_index_1] - y_pred_rf_contrast_back[i_temp_patient_index_1]

            # store i_temp_delta_prob in a data frame
            #       data frame row is patient and column is feature
            # pdb.set_trace()
            df_class_contrastive_ml_2comb[i_temp_patient_index_1, i_temp_row_counter] = i_temp_delta_prob
            
            # increment temporary row counter
            i_temp_row_counter = i_temp_row_counter + 1

            # restore test set
            x_test_orig = x_test_original_archive

            #  end for loop

    # restore test set
    x_test_orig = x_test_original_archive

    
    # save data frame of class-contrastive results to hard disk
    # this will be read from R
    np.savetxt('df_class_contrastive_ml_2comb_from_python.csv', df_class_contrastive_ml_2comb, delimiter = '\t')



    pdb.set_trace()


    ##################################################################
    # Automated machine learning using TPOT
    ################################################################## 	
    # reformat for TPOT (numpy arrays)
    # x_train_array = np.array(x_train)
    # y_train_array = np.array(y_train)

    # TODO: change this to all data
    #X_train, X_test, Y_train, Y_test = train_test_split(x_train_array_concat,
    #                                                    y_train_array, train_size=0.75, test_size=0.25)

    tpot = TPOTClassifier(generations=10, population_size=500, verbosity=2, n_jobs=2) #scoring='f1', n_jobs=2) scoring='roc_auc'	
    # 'accuracy', 'adjusted_rand_score', 'average_precision', 'balanced_accuracy', 'f1', 'f1_macro',
    #'f1_micro', 'f1_samples', 'f1_weighted', 'neg_log_loss','precision', 'precision_macro', 'precision_micro',
    #'precision_samples', 'precision_weighted', 'recall', 'recall_macro', 'recall_micro', 'recall_samples', 'recall_weighted', 'roc_auc'
	# See more at 
    # https://epistasislab.github.io/tpot/api/
 	
    tpot.fit(x_train, y_train)
    # TODO: use tpoit.predeict and then use roc_curve to calculate AUC and score should be accuracy	
    print("Score from a baseline automated machine learning algorithm: ", tpot.score(x_test, y_test))


    print("\n ***************************************** \n")
    print("   All tasks successfully completed \n")
    print(" ***************************************** \n")


    # TODO: custom loss functions
    # MOST IMPORTANT https://keras.io/examples/variational_autoencoder/
    #
    # https://www.tensorflow.org/api_docs/python/tf/keras/losses/KLD
    # https://www.tensorflow.org/api_docs/python/tf/keras/losses/binary_crossentropy
    # https://github.com/tensorflow/tensorflow/blob/2b96f3662bd776e277f86997659e61046b56c315/tensorflow/python/keras/losses.py#L1492
    # https://towardsdatascience.com/how-to-create-a-custom-loss-function-keras-3a89156ec69b
    # https://heartbeat.fritz.ai/how-to-create-a-custom-loss-function-in-keras-637bd312e9ab
    # https://www.google.com/search?client=firefox-b-d&q=call+custom+loss+function+keras
    # https://www.tensorflow.org/tutorials/generative/cvae
    # https://www.google.com/search?client=firefox-b-d&q=elbo+loss+keras
    # 
