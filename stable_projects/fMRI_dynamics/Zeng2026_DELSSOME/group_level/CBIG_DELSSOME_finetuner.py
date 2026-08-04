"""
Written by Tianchu Zeng and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import argparse
import configparser
from datetime import datetime
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import torch
import torch.nn as nn
import torch.utils.data as Data
import os
import gc

from models.CBIG_DELSSOME_model import MultiModalCostPredictor
from utils.CBIG_DELSSOME_pFIC_utils import CBIG_corr
from process.CBIG_DELSSOME_data_process import data_process_predictor
from process.CBIG_DELSSOME_datasets import DataSetPredictor


class PredictorTrainer:

    def __init__(self, config, finetuning=True, freeze_transformer=False, seed=None):

        self.seed = seed
        if seed is None:
            self.seed = np.random.randint(0, 1000000000)
        torch.manual_seed(self.seed)

        torch.set_default_dtype(torch.float64)
        if torch.cuda.is_available():
            self.device = 'cuda'
            torch.cuda.manual_seed(self.seed)
            print(torch.cuda.get_device_name())
        else:
            self.device = 'cpu'

        # Read config file
        data_config = config['data preparation']

        self.group_mats_path = data_config['group_mats_path']
        # The path to the group-level myelin, gradient, SC, FC and FCD-CDF matrices
        self.param_path = data_config['param_path']
        # The path to the parameter file with corresponding loss

        self.n_training_group = int(data_config['n_training_group'])
        self.n_val_group = int(data_config['n_val_group'])

        self.data_epochs = int(data_config.get('data_epochs', '100'))

        # Load data by batch, just for saving the memory
        self.train_group_list = np.arange(0, self.n_training_group)
        self.val_group_list = np.arange(0, self.n_val_group)

        predictor_training_config = config['predictor training']
        self.pre_epoch = 0
        # The number of previously trained epochs, for continuing training
        self.epochs = int(predictor_training_config['epochs'])
        # The number of epochs for training
        self.batch_size = int(predictor_training_config['batch_size'])
        # The batch size for training

        self.lr = float(predictor_training_config['lr'])
        # The learning rate for training
        self.gamma = float(predictor_training_config['gamma'])
        # The gamma for the learning rate scheduler

        model_save_dir = predictor_training_config['predictor_save_dir']
        os.makedirs(model_save_dir, exist_ok=True)
        self.model_save_path = os.path.join(model_save_dir, 'predictor_model.pth')

        predictor_model_config = config['predictor model']
        n_regions = int(predictor_model_config['n_regions'])
        d_param = int(predictor_model_config['d_param'])
        d = int(predictor_model_config['d'])

        self.model = MultiModalCostPredictor(n_regions=n_regions, d_param=d_param, d=d)

        if self.pre_epoch > 0:
            prev_model = torch.load(self.model_save_path)
            self.model.load_state_dict(prev_model['model_state_dict'])
            print("Loaded previous model state dict.")
        elif finetuning:
            # Fine-tuning starts from a pre-trained predictor. When
            # ``freeze_transformer`` is False all parameters are fine-tuned (pFIC
            # setup); when True the transformer layers are frozen and only the
            # remaining parameters are updated (pMFM setup). Useful when adapting
            # to a new subject group with limited training data.
            print("Finetuning the model with pre-trained weights...")
            pre_trained_model_path = predictor_training_config['pre_trained_model_path']
            if not os.path.exists(pre_trained_model_path):
                raise FileNotFoundError(f"Pre-trained model not found at {pre_trained_model_path}.")
            self.model.load_state_dict(
                torch.load(pre_trained_model_path, map_location=torch.device(self.device))['model_state_dict'])
            print(f"Loaded pre-trained model state dict from {pre_trained_model_path}.")

            if freeze_transformer:
                # Freeze the transformer layers, fine-tune everything else.
                for name, param in self.model.named_parameters():
                    param.requires_grad = not name.startswith("transformer")
                print("Froze transformer layers for fine-tuning.")

        if torch.cuda.is_available():
            self.model = self.model.cuda()

        # ``filter`` is a no-op when nothing is frozen (all params trainable), so
        # this covers both the frozen (pMFM) and unfrozen (pFIC/Hopf) cases.
        self.optimizer = torch.optim.Adam(filter(lambda p: p.requires_grad, self.model.parameters()), lr=self.lr)
        self.scheduler = torch.optim.lr_scheduler.ExponentialLR(self.optimizer, gamma=self.gamma)
        self.loss_func = nn.MSELoss()

        # Load data
        train_data = data_process_predictor(mode='train',
                                            n_epoch=self.data_epochs,
                                            group_mats_path=self.group_mats_path,
                                            param_path=self.param_path,
                                            group_list=self.train_group_list)
        self.train_dataset = DataSetPredictor(*train_data)
        self.train_loader = Data.DataLoader(self.train_dataset,
                                            self.batch_size,
                                            shuffle=True,
                                            generator=torch.Generator(device='cpu'))

        val_data = data_process_predictor(mode='val',
                                          n_epoch=self.data_epochs,
                                          group_mats_path=self.group_mats_path,
                                          param_path=self.param_path,
                                          group_list=self.val_group_list)
        self.val_dataset = DataSetPredictor(*val_data)
        self.val_loader = Data.DataLoader(self.val_dataset,
                                          self.batch_size,
                                          shuffle=False,
                                          generator=torch.Generator(device='cpu'))

        print("DELSSOME FC+FCD cost predictor trainer has been successfully initialized.")

    def train(self, need_early_stop=False):
        """The main function for training the DELSSOME predictor.

        Args:
            need_early_stop (bool, optional): Whether to activate early stop. Defaults to False.

        Return 0 if no error occurs.

        """

        train_loss_list = []
        val_loss_list = []
        min_val_loss = 1000
        ascending_count = 0
        for epoch in range(self.pre_epoch, self.epochs):
            print('Epoch: [{}/{}], batch size: {}, lr: {}'.format(epoch + 1, self.epochs, self.batch_size,
                                                                  self.optimizer.param_groups[0]['lr']))
            start_time = datetime.now()
            train_loss_epoch = self._train_one_epoch_by_all(epoch)
            lapse_time = datetime.now() - start_time
            print(f'Average train loss: {train_loss_epoch}; Time cost: {lapse_time}')
            train_loss_list.append(train_loss_epoch)
            print('Start evaluating...')
            start_time = datetime.now()
            val_loss_epoch = self._evaluate_by_all()
            lapse_time = datetime.now() - start_time
            print(f'Val loss: {val_loss_epoch}; Time cost: {lapse_time}')
            val_loss_list.append(val_loss_epoch)
            if val_loss_epoch < min_val_loss:
                min_val_loss = val_loss_epoch
                print(f"Get a minimal validation loss ever. Start saving to {self.model_save_path}...")
                torch.save(
                    {
                        'model_state_dict': self.model.state_dict(),
                        'seed': self.seed,
                        'train_loss': train_loss_list,
                        'val_loss': val_loss_list
                    }, self.model_save_path)
                print("Saved successfully.")
            if need_early_stop:
                if epoch > self.pre_epoch + 5:
                    if val_loss_list[epoch - self.pre_epoch] > val_loss_list[epoch - self.pre_epoch - 1]:
                        ascending_count += 1
                    else:
                        ascending_count = 0
                    if ascending_count >= 3:
                        print("The loss is ascending. Early break.")
                        break
        return 0

    def _evaluate_by_all(self):
        self.model.eval()
        loss_record = 0
        with torch.no_grad():

            torch.cuda.empty_cache()
            gc.collect()

            loss_total = 0

            for parameter, ground_loss, sc_mat, emp_fc, emp_fcd in self.val_loader:

                if torch.cuda.is_available():
                    parameter = parameter.cuda()
                    sc_mat = sc_mat.cuda()
                    emp_fc = emp_fc.cuda()
                    emp_fcd = emp_fcd.cuda()
                    ground_loss = ground_loss.cuda()

                pred_loss = self.model(parameter, sc_mat, emp_fc, emp_fcd)

                loss = self.loss_func(pred_loss, ground_loss)
                loss_total += loss.item()

            loss_record += loss_total

            return loss_record / self.n_val_group

    def _train_one_epoch_by_all(self, epoch):
        self.model.train()
        loss_record = 0

        # Clear memory before loading new data
        torch.cuda.empty_cache()
        gc.collect()

        train_steps = len(self.train_loader)

        for i, (parameter, ground_loss, sc_mat, emp_fc, emp_fcd) in enumerate(self.train_loader):
            if parameter.shape[
                    0] == 1:  # To avoid batch size equals to 1 because we use BatchNorm in deep learning models
                continue

            if torch.cuda.is_available():
                parameter = parameter.cuda()
                sc_mat = sc_mat.cuda()
                emp_fc = emp_fc.cuda()
                emp_fcd = emp_fcd.cuda()
                ground_loss = ground_loss.cuda()

            self.optimizer.zero_grad()
            pred_loss = self.model(parameter, sc_mat, emp_fc, emp_fcd)

            loss = self.loss_func(pred_loss, ground_loss)
            loss.backward()
            self.optimizer.step()
            loss_record += loss.item()
            if (i + 1) % 100 == 0:
                print("Epoch [{}/{}], step [{}/{}], loss:{:.4f}".format(epoch + 1, self.epochs, i + 1, train_steps,
                                                                        loss.item()))

        self.scheduler.step()
        return loss_record / self.n_training_group


class PredictorTester:

    def __init__(self, config):

        torch.set_default_dtype(torch.float64)
        if torch.cuda.is_available():
            self.device = 'cuda'
        else:
            self.device = 'cpu'

        # Read config file
        data_config = config['data preparation']
        self.group_mats_path = data_config['group_mats_path']
        # The path to the group-level myelin, gradient, SC, FC and FCD-CDF matrices
        self.param_path = data_config['param_path']
        # The path to the parameter file with corresponding loss
        self.n_test_group = int(data_config['n_test_group'])

        self.data_epochs = int(data_config.get('data_epochs', '100'))
        self.batch_size = int(config['predictor training']['batch_size'])

        self.test_group_list = np.arange(0, self.n_test_group)
        self.test_group_batch = 1

        model_save_dir = config['predictor training']['predictor_save_dir']
        self.model_save_path = os.path.join(model_save_dir, 'predictor_model.pth')
        # Save record file path
        self.record_save_dir = os.path.join(model_save_dir, 'predictor_record')
        os.makedirs(self.record_save_dir, exist_ok=True)
        self.record_path = os.path.join(self.record_save_dir, 'predictor_record.pth')
        # Save figure path
        self.figure_save_dir = os.path.join(self.record_save_dir, 'figures')

        predictor_model_config = config['predictor model']
        n_regions = int(predictor_model_config['n_regions'])
        d_param = int(predictor_model_config['d_param'])
        d = int(predictor_model_config['d'])

        self.model = MultiModalCostPredictor(n_regions=n_regions, d_param=d_param, d=d)

        self.model.load_state_dict(
            torch.load(self.model_save_path, map_location=torch.device(self.device))['model_state_dict'])
        if torch.cuda.is_available():
            self.model = self.model.cuda()
        print(f'Model successfully loaded from {self.model_save_path}.')

        print("DELSSOME FC+FCD cost predictor tester has been successfully initialized.")

    def test(self, need_boxplot=False, need_loss_figure=False):
        """The main function for testing the DELSSOME predictor.

        Args:
            need_boxplot (bool, optional): Whether to plot the boxplot. Defaults to False.
            need_loss_figure (bool, optional): Whether to plot the loss figure. Defaults to False.

        Return 0 if no error occurs.

        """

        self.need_boxplot = need_boxplot
        self.need_loss_figure = need_loss_figure

        print(" -- Start testing -- ")
        if need_boxplot or need_loss_figure:
            if not os.path.exists(self.figure_save_dir):
                os.mkdir(self.figure_save_dir)
            print(f"The figures will be saved to {self.figure_save_dir}.")
        self._test_by_set()
        print(" -- Test Done -- ")
        return 0

    def _test_by_set(self):
        self.model.eval()
        with torch.no_grad():
            self.mse_mean = torch.zeros(self.n_test_group, 3)
            self.mse_std = torch.zeros(self.n_test_group, 3)
            self.corr_pred_ground = torch.zeros(self.n_test_group, 3)
            for group_nbr in range(self.n_test_group):
                print(f"Set [{group_nbr + 1}/{self.n_test_group}]")
                self._test_one_set(group_nbr)
            # End for all sets

            # Results
            print("Start saving data...")
            print(f"Record path: {self.record_path}")
            torch.save({
                'mse_mean': self.mse_mean,
                'mse_std': self.mse_std,
                'corr_pred_ground': self.corr_pred_ground
            }, self.record_path)
            print("Successfully saved the record file.")

            mean_mse_mean = torch.mean(self.mse_mean, dim=0)
            mean_mse_std = torch.mean(self.mse_std, dim=0)
            mean_corr_pred_ground = torch.mean(self.corr_pred_ground, dim=0)
            print(f"The MSE loss mean in whole test set is {mean_mse_mean.numpy()}")
            print(f"The MSE loss std in whole test set is {mean_mse_std.numpy()}")
            print(
                f"The correlation between predicted loss and ground loss in test set is {mean_corr_pred_ground.numpy()}"
            )
            return 0

    def _test_one_set(self, group_nbr):
        self.model.eval()
        with torch.no_grad():

            # Clear memory before loading new data
            torch.cuda.empty_cache()
            gc.collect()

            group_batch_list = self.test_group_list[group_nbr:group_nbr + self.test_group_batch]
            test_data = data_process_predictor(mode='test',
                                               n_epoch=self.data_epochs,
                                               group_mats_path=self.group_mats_path,
                                               param_path=self.param_path,
                                               group_list=group_batch_list)
            test_dataset = DataSetPredictor(*test_data)

            # Put all data into one batch, so set the batch size to the number of parameter sets
            test_loader = Data.DataLoader(dataset=test_dataset,
                                          batch_size=test_dataset.parameter_sets.shape[1],
                                          shuffle=False,
                                          generator=torch.Generator(device='cpu'))
            for parameter, ground_loss, sc_mat, emp_fc, emp_fcd in test_loader:

                if torch.cuda.is_available():
                    parameter = parameter.cuda()
                    sc_mat = sc_mat.cuda()
                    emp_fc = emp_fc.cuda()
                    emp_fcd = emp_fcd.cuda()
                    ground_loss = ground_loss.cuda()

                pred_loss = self.model(parameter, sc_mat, emp_fc, emp_fcd)

            if torch.cuda.is_available():
                pred_loss = pred_loss.cpu()
                ground_loss = ground_loss.cpu()
            mse_cur_set = torch.square(pred_loss - ground_loss)
            self.mse_mean[group_nbr] = torch.mean(mse_cur_set, dim=0)
            self.mse_std[group_nbr] = torch.std(mse_cur_set, dim=0)
            for j in range(3):
                self.corr_pred_ground[group_nbr, j] = CBIG_corr(pred_loss[:, j].unsqueeze(1),
                                                                ground_loss[:, j].unsqueeze(1))
            print(f"The MSE loss mean in test set {group_nbr} is {self.mse_mean[group_nbr].numpy()}")
            print(f"The MSE loss std in test set {group_nbr} is {self.mse_std[group_nbr].numpy()}")
            print((f"The correlation between predicted loss and ground loss in test set {group_nbr}"
                   f"is {self.corr_pred_ground[group_nbr].numpy()}"))

            # Plots
            if self.need_boxplot:
                print('Drawing box plot...')
                plt.figure()
                plt.boxplot(mse_cur_set.numpy(), labels=['Corr_loss', 'L1_loss', 'KS_loss'], showfliers=False)
                plt.xlabel('Loss type')
                plt.ylabel('MSE loss')
                figure_title = f'boxplot_set{group_nbr}'
                plt.title(figure_title)
                plt.savefig(os.path.join(self.figure_save_dir, f'{figure_title}.svg'))
                plt.close()
                print('Figure saved.')
            if self.need_loss_figure:
                print('Drawing loss figure...')
                torch.save({
                    'pred_loss': pred_loss,
                    'ground_loss': ground_loss
                }, os.path.join(self.figure_save_dir, f'loss_set{group_nbr}.pth'))
                print('Data saved.')

    def plot_loss_figure_one_set(self, pred_loss, ground_loss, set_nbr):
        figure_name = ['corr_loss', 'l1_loss', 'ks_loss']
        for i in range(len(figure_name)):
            r = CBIG_corr(pred_loss[:, i].unsqueeze(1), ground_loss[:, i].unsqueeze(1))
            plt.figure()
            sns.regplot(
                x=pred_loss[:, i].numpy(),
                y=ground_loss[:, i].numpy(),
                # yapf: disable
                scatter_kws={
                    'color': '#696969',
                    's': 10
                },
                # yapf: enable
                line_kws={'color': 'red'},
                order=1)
            plt.xlabel("Predicted loss")
            plt.ylabel("Ground loss")
            figure_title = f'{figure_name[i]}_set{set_nbr}'
            plt.title(f'{figure_title}, r={r.item():.4f}')
            plt.savefig(os.path.join(self.figure_save_dir, f'{figure_title}.svg'))
            plt.close()
            print(f"Figure {i + 1} saved.")
        return 0


def train_main(config_path, seed=None):
    # Read config file
    config = configparser.ConfigParser()
    config.read(config_path)

    predictor_training_config = config['predictor training']
    # `finetuning`: load a pre-trained predictor before training (pMFM) vs. train
    # from scratch (pFIC/Hopf, which the authors found works better).
    # `freeze_transformer`: freeze the transformer layers while fine-tuning (pMFM).
    # Both default to False to preserve the original pFIC behavior.
    finetuning = predictor_training_config.getboolean('finetuning', fallback=False)
    freeze_transformer = predictor_training_config.getboolean('freeze_transformer', fallback=False)

    trainer = PredictorTrainer(config, finetuning=finetuning, freeze_transformer=freeze_transformer, seed=seed)
    trainer.train()

    return 0


def test_main(config_path):
    # Read config file
    config = configparser.ConfigParser()
    config.read(config_path)

    tester = PredictorTester(config)
    tester.test(need_boxplot=False, need_loss_figure=True)

    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='DELSSOME predictor finetuner (pFIC / pMFM / Hopf)')
    parser.add_argument('--config_path',
                        type=str,
                        default='../examples/config/example_finetune_pFIC.ini',
                        help='The path to config file')
    parser.add_argument('--seed', type=int, default=None, help='The random seed for training')
    args = parser.parse_args()
    train_main(args.config_path, seed=args.seed)
    test_main(args.config_path)
