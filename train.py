import tensorflow as tf
from tensorflow.python.ops import math_ops
from tensorflow.python import ops
import keras
import numpy as np
import tables
from keras.callbacks import ModelCheckpoint, EarlyStopping
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import argparse
from models import dense
import math
import keras.backend as K
import uproot

def huber_loss(y_true, y_pred, delta=1.0):
    error = y_pred - y_true
    abs_error = K.abs(error)
    quadratic = K.minimum(abs_error, delta)
    linear = abs_error - quadratic
    return 0.5 * K.square(quadratic) + delta * linear

def mean_absolute_relative_error(y_true, y_pred):
    if not K.is_tensor(y_pred):
        y_pred = K.constant(y_pred)
    y_true = K.cast(y_true, y_pred.dtype)
    diff = K.abs((y_true - y_pred) / K.clip(K.abs(y_true),
                                            K.epsilon(),
                                            None))
    return K.mean(diff, axis=-1)

def mean_squared_relative_error(y_true, y_pred):
    if not K.is_tensor(y_pred):
        y_pred = K.constant(y_pred)
    y_true = K.cast(y_true, y_pred.dtype)
    diff = K.square((y_true - y_pred) / K.clip(K.square(y_true),
                                            K.epsilon(),
                                            None))
    return K.mean(diff, axis=-1)

def get_features_targets(file_name, features, targets, spectators):
    # load file
    h5file = tables.open_file(file_name, "r")
    nevents = getattr(h5file.root,features[0]).shape[0]
    ntargets = len(targets)
    nfeatures = len(features)
    nspec = len(spectators)

    # allocate arrays
    feature_array = np.zeros((nevents,nfeatures))
    target_array = np.zeros((nevents,ntargets))
    spec_array = np.zeros((nevents,nspec))

    # load feature arrays
    for (i, feat) in enumerate(features):
        feature_array[:,i] = getattr(h5file.root,feat)[:]
    # load target arrays
    for (i, targ) in enumerate(targets):
        target_array[:,i] = getattr(h5file.root,targ)[:]
    for (i, targ) in enumerate(spectators):
        spec_array[:,i] = getattr(h5file.root,targ)[:]

    h5file.close()
    return feature_array,target_array,spec_array

def write_outputs(file_name, features, targets, spectators):
    # write
    h5file = tables.open_file(file_name, "r")
    nevents = getattr(h5file.root,features[0]).shape[0]
    ntargets = len(targets)
    nfeatures = len(features)
    nspec = len(spectators)

    # allocate arrays
    feature_array = np.zeros((nevents,nfeatures))
    target_array = np.zeros((nevents,ntargets))
    spec_array = np.zeros((nevents,nspec))

    # load feature arrays
    for (i, feat) in enumerate(features):
        feature_array[:,i] = getattr(h5file.root,feat)[:]
    # load target arrays
    for (i, targ) in enumerate(targets):
        target_array[:,i] = getattr(h5file.root,targ)[:]
    for (i, targ) in enumerate(spectators):
        spec_array[:,i] = getattr(h5file.root,targ)[:]

    h5file.close()

def main(args):
    file_path = "samples/input_MET_DY.h5"

    features = [
        "L1CaloMet_para_puppi","L1CaloMet_perp_puppi",
        "L1PFMet_para_puppi","L1PFMet_perp_puppi",
        "L1PuppiMet_para_puppi","L1PuppiMet_perp_puppi",
        "CaloHTMiss_para_puppi","CaloHTMiss_perp_puppi",
        "PFHTMiss_para_puppi","PFHTMiss_perp_puppi",
        "PuppiHTMiss_para_puppi","PuppiHTMiss_perp_puppi",
    ]

    targets = ["genMet_para_puppi","genMet_perp_puppi",]
    spectators = ["event","zPt","zPhi","L1PuppiMet","L1PuppiMetPhi"]

    feature_array, target_array, spectator_array = get_features_targets(file_path, features, targets, spectators)
    nevents = feature_array.shape[0]
    nfeatures = feature_array.shape[1]
    ntargets = target_array.shape[1]
    nspec = spectator_array.shape[1]
    print("Loaded {} events with {} input features, {} targets, and {} spectator variables".format(
        nevents,nfeatures,ntargets,nspec))

    # fit keras model
    X = feature_array
    y = target_array

    fulllen = nevents
    tv_frac = 0.10
    tv_num = math.ceil(fulllen*tv_frac)
    splits = np.cumsum([fulllen-2*tv_num,tv_num,tv_num])
    splits = [int(s) for s in splits]

    X_train = X[0:splits[0]]
    X_val = X[splits[1]:splits[2]]
    X_test = X[splits[0]:splits[1]]

    y_train = y[0:splits[0]]
    y_val = y[splits[1]:splits[2]]
    y_test = y[splits[0]:splits[1]]

    spec_train = spectator_array[0:splits[0]]
    spec_val = spectator_array[splits[1]:splits[2]]
    spec_test = spectator_array[splits[0]:splits[1]]

    def loss_para(y_true, y_pred):
        return K.mean(math_ops.square(y_pred - y_true), axis=-1)
    def loss_perp(y_true, y_pred):
        return K.mean(math_ops.square(y_pred - y_true), axis=-1)

    keras_model = dense(nfeatures, ntargets)

    keras_model.compile(optimizer="adam", loss=[loss_para, loss_perp], 
                        loss_weights = [1., 1.], metrics=["mean_absolute_error"])
    print(keras_model.summary())

    early_stopping = EarlyStopping(monitor="val_loss", patience=10)
    model_checkpoint = ModelCheckpoint("keras_model_best.h5", monitor="val_loss", save_best_only=True)
    callbacks = [early_stopping, model_checkpoint]

    keras_model.fit(X_train, [y_train[:,:1], y_train[:,1:]], batch_size=1024, 
                    epochs=100, validation_data=(X_val, [y_val[:,:1], y_val[:,1:]]), shuffle=True,
                    callbacks = callbacks)

    keras_model.load_weights("keras_model_best.h5")
    
    predict_test = keras_model.predict(X_test)
    predict_test = np.concatenate(predict_test,axis=1)
    # print(y_test)
    # print(predict_test)

    def print_res(gen_met, predict_met, name="Met_res.pdf"):
        plt.figure()          
        #rel_err = (predict_met - gen_met)/np.clip(gen_met, 1e-6, None)
        #plt.hist(rel_err, bins=np.linspace(-1., 1., 50+1))
        #plt.hist(err, bins=np.linspace(-1., 1., 50+1))
        err = (predict_met - gen_met)
        plt.hist(err, bins=np.linspace(-100., 100., 50+1))
        plt.xlabel("Rel. err.")
        plt.ylabel("Events")
        plt.figtext(0.25, 0.90,"CMS",fontweight="bold", wrap=True, horizontalalignment="right", fontsize=14)
        plt.figtext(0.35, 0.90,"preliminary", style="italic", wrap=True, horizontalalignment="center", fontsize=14) 
        plt.savefig(name)
	
    print_res(y_test[:,0], predict_test[:,0], name = "MET_para_res.pdf")
    print_res(y_test[:,1], predict_test[:,1], name = "MET_perp_res.pdf")

    #generate other outputs
    predict_val = keras_model.predict(X_val)
    predict_val = np.concatenate(predict_val,axis=1)
    predict_train = keras_model.predict(X_train)
    predict_train = np.concatenate(predict_train,axis=1)

    # output prediction, event, type (test/train/valid)
    out_pred = np.concatenate((predict_train,predict_val,predict_test))
    out_type = np.concatenate((np.zeros(predict_train.shape[0]),
                               np.ones(predict_val.shape[0]),
                               2*np.ones(predict_test.shape[0])) )
    out_spec = np.concatenate((spec_train,spec_val,spec_test))
    #print(out_spec.shape)
    # print(spec_train.shape)
    # print(spec_val.shape)
    # print(spec_test.shape)
    #spectators = ["event","zPt","zPhi","L1PuppiMet","L1PuppiMetPhi"]


    # def _write_carray(a, h5file, name, group_path="/", **kwargs):
    #     h5file.create_carray(group_path, name, obj=a, filters=filters, createparents=True, **kwargs)
    # _write_carray(df_met[k].values, h5file, name=k.replace("[","").replace("]",""))

    # with tables.open_file("ouput_MET_DY.h5", mode="w") as hf:
    #     filters = tables.Filters(complevel=7, complib="blosc")
    #     hf.create_carray("/","pred_para",out_pred[:,0],filters=filters)
    #     hf.create_carray("/","pred_perp", out_pred[:,1],filters=filters)
    #     hf.create_carray("/","event", out_evt,filters=filters)
    #     hf.create_carray("/","type", out_type,filters=filters)

    import h5py
    hf = h5py.File("samples/output_MET_DY.h5", "w")
    hf.create_dataset("pred_para_puppi", data = out_pred[:,0])
    hf.create_dataset("pred_perp_puppi", data = out_pred[:,1])
    #hf.create_dataset("event",     data = out_evt)
    hf.create_dataset("type",      data = out_type)
    for isp, sp in enumerate(spectators):
        hf.create_dataset(sp, data = out_spec[:,isp])

    # can also save as a root file from h5
    f = uproot.recreate("samples/output_MET_DY.root")
    f["Events"] = uproot.newtree( hf )
    f["Events"].extend( {k:np.array(hf[k]) for k in hf.keys()} )
    f.close()
        
    hf.close()

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
        
    args = parser.parse_args()
    main(args)
