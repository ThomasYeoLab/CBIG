from __future__ import print_function
import time
import pickle

import cbig.Nguyen2020.train as train
import cbig.Nguyen2020.predict as predict
import cbig.Nguyen2020.evaluation as ev
import cbig.Nguyen2020.misc as misc


if __name__ == '__main__':
    timestamp = time.strftime('.%m%d-%H%M%S', time.localtime())

    # Parse command line arguments
    train_args, _ = train.get_arg_parser().parse_known_args()
    val_args, _ = predict.get_arg_parser().parse_known_args()
    ev_args, _ = ev.get_arg_parser().parse_known_args()

    # Rename to avoid writing to the same file during HORD optimization
    checkpoint = train_args.checkpoint + timestamp + '.pt'
    setattr(train_args, 'checkpoint', checkpoint)
    setattr(val_args, 'checkpoint', checkpoint)
    prediction = val_args.prediction + timestamp + '.csv'
    setattr(val_args, 'prediction', prediction)
    setattr(ev_args, 'prediction', prediction)

    train.train(train_args)
    predict.main(val_args)
    result = ev.eval_submission(misc.read_csv(ev_args.reference), misc.read_csv(ev_args.prediction))

    with open(val_args.data, 'rb') as fhandler:
        data = pickle.load(fhandler)

    # Normalize ADAS and Vent result for HORD optimization
    adas = result['adasMAE'] / data['stds']['ADAS13']
    vent = result['ventsMAE'] / data['VentICVstd']

    print(prediction, result['mAUC'], result['bca'], adas, vent, sep=',')
