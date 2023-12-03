from TransgenicMouseData import TransgenicMouseData

import numpy as np
import pandas as pd
import openpyxl as op
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.optimize import curve_fit
from datetime import date
import math
import json 
import os

cwd = "/Users/veramazeeva/Library/CloudStorage/OneDrive-Personal/Carnegie Mellon/CompBio512/project"

class TumorGrowth():

    def __init__(self, untreated_only=True):

        self.id2iterator = {0: 'fe', 1: 'lf', 2: 'be'}
        self.iterator2id = {'fe': 0, 'lf': 1, 'be': 2}

        mice = TransgenicMouseData()

        self.mouse2data = dict((m, mice.mouse2data[m]) for m in mice.untreated_mice)

        if not untreated_only:
            self.mouse2data = mice.mouse2data
        
        for mouse_id in self.mouse2data.keys():
            self.mouse2data[mouse_id]['Day'] = [0]
            for i in range(1, len(self.mouse2data[mouse_id]['Date'])):
                self.mouse2data[mouse_id]['Day'] += [(self.mouse2data[mouse_id]['Date'][i] - self.mouse2data[mouse_id]['Date'][i-1]).days + self.mouse2data[mouse_id]['Day'][-1]]

        self.mouse2model = dict()
        self.mouse2fit = dict()

        for mouse_id, mouse_data in self.mouse2data.items():
            self.mouse2model[mouse_id] = dict()
            self.mouse2fit[mouse_id] = dict()
            self.__init__models(mouse_id, mouse_data['Volume'])
    
    def __init__models(self, mouse, volume_data):
        self.mouse2model[mouse]['Logistic'], self.mouse2fit[mouse]['Logistic'] = self.fit_logistic(volume_data)
        self.mouse2model[mouse]['Gompertz'], self.mouse2fit[mouse]['Gompertz'] = self.fit_gompertz(volume_data)
        self.mouse2model[mouse]['Bertalanffy'], self.mouse2fit[mouse]['Bertalanffy'] = self.fit_bertalanffy(volume_data)
    

    ################################ Logistic ################################

    def fit_logistic(self, data):
        N = len(data)
         
        def logistic_growth(V, a, b):
            return a * V - b * V**2
         
        def logistic_fe(T, a, b):
            n = len(T)
            solution = np.zeros_like(T)
            solution[0] = data[0]
            for i in range(1, n):
                solution[i] = solution[i-1] + logistic_growth(solution[i-1], a, b)
            return solution 
        
        def logistic_lf(T, a, b):
            n = len(T)
            solution = np.zeros_like(T)
            solution[0] = data[0]
            solution[1] = solution[0] + logistic_growth(solution[0], a, b)
            for i in range(2, n):
                solution[i] = solution[i-2] + logistic_growth(solution[i-1], a, b)
            return solution 
        
        def f(x, x_prev, a, b):
             return (x_prev - logistic_growth(x, a, b))

        def logistic_be(T, a, b):
            n = len(T)
            solution = np.zeros_like(T)
            solution[0] = data[0]
            for i in range(1, n):
                prev = solution[i-1]
                solution[i], = fsolve(f, prev, args=(a, b, prev))
                solution[i] = max(solution[i], 0)
            return solution 

        guess = [0, 0]
        bounds = ([0, 0,], [np.inf, np.inf])

        T = [i for i in range(N)]

        params_fe, _ = curve_fit(logistic_fe, T, data, p0=guess, bounds=bounds, maxfev=20000)
        params_lf, _ = curve_fit(logistic_lf, T, data, p0=guess, bounds=bounds, maxfev=20000)
        params_be, _ = curve_fit(logistic_be, T, data, p0=guess, bounds=bounds, maxfev=20000)
        
        def logistic_solution(V0, params, T):
            (a, b) = params
            solution = np.zeros_like(T)
            solution[0] = V0
            for i in range(1, len(T)):
                s = (V0 * (a/b) * np.exp(a*i)) / (a/b - V0 + V0 * np.exp(a*i))
                solution[i] = s
            return solution 
        
        fitted_fe = logistic_solution(data[0], params_fe, T)
        fitted_lf = logistic_solution(data[0], params_lf, T)
        fitted_be = logistic_solution(data[0], params_be, T)

        two = [2] * N
        nmse_fe = sum(np.power(fitted_fe - data[:len(fitted_fe)], two)) / sum(np.power(data[:len(fitted_fe)], two))
        nmse_lf = sum(np.power(fitted_lf - data[:len(fitted_lf)], two)) / sum(np.power(data[:len(fitted_lf)], two))
        nmse_be = sum(np.power(fitted_be - data[:len(fitted_be)], two)) / sum(np.power(data[:len(fitted_be)], two))
        
        Model = {'params': [params_fe, params_lf, params_be], 'nmse': [nmse_fe, nmse_lf, nmse_be]}
        Fits = [fitted_fe, fitted_lf, fitted_be]

        return Model, Fits

    ################################ Gompertz ################################

    def fit_gompertz(self, data):
        N = len(data)
         
        def gompertz_growth(V, a, b):
            return a * V - b * V * np.log(V)
         
        def gompertz_fe(T, a, b):
            n = len(T)
            solution = np.zeros_like(T)
            solution[0] = data[0]
            for i in range(1, n):
                solution[i] = solution[i-1] + gompertz_growth(solution[i-1], a, b)
            return solution 
        
        def gompertz_lf(T, a, b):
            n = len(T)
            solution = np.zeros_like(T)
            solution[0] = data[0]
            solution[1] = solution[0] + gompertz_growth(solution[0], a, b)
            for i in range(2, n):
                solution[i] = solution[i-2] + gompertz_growth(solution[i-1], a, b)
            return solution 
        
        def f(x, x_prev, a, b):
             return (x_prev - gompertz_growth(x, a, b))

        def gompertz_be(T, a, b):
            n = len(T)
            solution = np.zeros_like(T)
            solution[0] = data[0]
            for i in range(1, n):
                prev = solution[i-1]
                solution[i], = fsolve(f, prev, args=(a, b, prev))
                solution[i] = max(solution[i], 0)
            return solution 

        guess = [0, 0]
        bounds = ([0, 0,], [np.inf, np.inf])

        T = [i for i in range(N)]

        params_fe, _ = curve_fit(gompertz_fe, T, data, p0=guess, bounds=bounds, maxfev=20000)
        params_lf, _ = curve_fit(gompertz_lf, T, data, p0=guess, bounds=bounds, maxfev=20000)
        params_be, _ = curve_fit(gompertz_be, T, data, p0=guess, bounds=bounds, maxfev=20000)

        def gompertz_solution(V0, params, T):
            (a, b) = params
            solution = np.zeros_like(T)
            solution[0] = V0
            for i in range(1, len(T)):
                s = np.exp((np.log(V0) - a/b)*np.exp(-b*i)+a/b)
                solution[i] = s
            return solution 
        
        fitted_fe = gompertz_solution(data[0], params_fe, T)
        fitted_lf = gompertz_solution(data[0], params_lf, T)
        fitted_be = gompertz_solution(data[0], params_be, T)

        two = [2] * N
        nmse_fe = sum(np.power(fitted_fe - data[:len(fitted_fe)], two)) / sum(np.power(data[:len(fitted_fe)], two))
        nmse_lf = sum(np.power(fitted_lf - data[:len(fitted_lf)], two)) / sum(np.power(data[:len(fitted_lf)], two))
        nmse_be = sum(np.power(fitted_be - data[:len(fitted_be)], two)) / sum(np.power(data[:len(fitted_be)], two))
            
        Model = {'params': [params_fe, params_lf, params_be], 'nmse': [nmse_fe, nmse_lf, nmse_be]}
        Fits = [fitted_fe, fitted_lf, fitted_be]

        return Model, Fits


     ################################ Bertalanffy ################################

    def fit_bertalanffy(self, data):
        N = len(data)
         
        def bertalanffy_growth(V, a, b):
            return a * V**(2/3) - b * V
         
        def bertalanffy_fe(T, a, b):
            n = len(T)
            solution = np.zeros_like(T)
            solution[0] = data[0]
            for i in range(1, n):
                solution[i] = solution[i-1] + bertalanffy_growth(solution[i-1], a, b)
            return solution 
        
        def bertalanffy_lf(T, a, b):
            n = len(T)
            solution = np.zeros_like(T)
            solution[0] = data[0]
            solution[1] = solution[0] + bertalanffy_growth(solution[0], a, b)
            for i in range(2, n):
                solution[i] = solution[i-2] + bertalanffy_growth(solution[i-1], a, b)
            return solution 
        
        def f(x, x_prev, a, b):
             return (x_prev - bertalanffy_growth(x, a, b))

        def bertalanffy_be(T, a, b):
            n = len(T)
            solution = np.zeros_like(T)
            solution[0] = data[0]
            for i in range(1, n):
                prev = solution[i-1]
                solution[i], = fsolve(f, prev, args=(a, b, prev))
                solution[i] = max(solution[i], 0)
            return solution 

        guess = [0, 0]
        bounds = ([0, 0,], [np.inf, np.inf])

        T = [i for i in range(N)]

        params_fe, _ = curve_fit(bertalanffy_fe, T, data, p0=guess, bounds=bounds, maxfev=20000)
        params_lf, _ = curve_fit(bertalanffy_lf, T, data, p0=guess, bounds=bounds, maxfev=20000)
        params_be, _ = curve_fit(bertalanffy_be, T, data, p0=guess, bounds=bounds, maxfev=20000)

        
        def bertalanffy_solution(V0, params, T):
            (a, b) = params
            solution = np.zeros_like(T)
            solution[0] = V0
            for i in range(1, len(T)):
                s = (a/b + (V0**(1/3) - a/b) * np.exp(-b*i/3))**3
                solution[i] = s
            return solution 
        
        fitted_fe = bertalanffy_solution(data[0], params_fe, T)
        fitted_lf = bertalanffy_solution(data[0], params_lf, T)
        fitted_be = bertalanffy_solution(data[0], params_be, T)

        two = [2] * N
        nmse_fe = sum(np.power(fitted_fe - data[:len(fitted_fe)], two)) / sum(np.power(data[:len(fitted_fe)], two))
        nmse_lf = sum(np.power(fitted_lf - data[:len(fitted_lf)], two)) / sum(np.power(data[:len(fitted_lf)], two))
        nmse_be = sum(np.power(fitted_be - data[:len(fitted_be)], two)) / sum(np.power(data[:len(fitted_be)], two))
            
        Model = {'params': [params_fe, params_lf, params_be], 'nmse': [nmse_fe, nmse_lf, nmse_be]}
        Fits = [fitted_fe, fitted_lf, fitted_be]

        return Model, Fits


    ################### Plotting ########################

    def plot(self, mouse_id, logistic_it, gompertz_it, bertalanffy_it):
        volume_data = self.mouse2data[mouse_id]['Volume']

        logistic_fit = self.mouse2fit[mouse_id]['Logistic'][self.iterator2id[logistic_it]]
        gompertz_fit = self.mouse2fit[mouse_id]['Gompertz'][self.iterator2id[gompertz_it]]
        bertalanffy_fit = self.mouse2fit[mouse_id]['Bertalanffy'][self.iterator2id[bertalanffy_it]]
        plt.figure()
        plt.plot([i for i in range(len(logistic_fit))], logistic_fit, label="logistic fitted curve")
        plt.plot([i for i in range(len(gompertz_fit))], gompertz_fit, label="gompertz fitted curve")
        plt.plot([i for i in range(len(bertalanffy_fit))], bertalanffy_fit, label="bertalanffy fitted curve")
        plt.plot([i for i in range(len(volume_data))], volume_data, label="tumor volume data")
        plt.legend()
        plt.title("Tumor Growth for Mouse " + mouse_id)
        plt.xlabel("Time")
        plt.ylabel("Volume")

        path = os.path.join(cwd, mouse_id + '.png')
        
        plt.savefig(path)

        #plt.show()






##################################### End of Class Definition ################################

def convert_all_models(M):
    new_model_dict = dict()
    for mouse_id, mouse_models in M.items():
        models = dict()
        for model in ['Logistic', 'Gompertz', 'Bertalanffy']:
            models[model] = dict()
            models[model]['params'] = []
            models[model]['nmse'] = []
            for p in mouse_models[model]['params']:
                models[model]['params'] += [tuple(p)]
            for e in mouse_models[model]['nmse']:
                models[model]['nmse'] += [e] 
        new_model_dict[mouse_id] = models
    return new_model_dict

def get_best_models(M):
    id2iterator = {0: 'fe', 1: 'lf', 2: 'be'}
    best_params = dict()
    best_models = dict()
    for mouse_id, mouse_models in M.items():
        best_params[mouse_id] = dict()
        best_model = ''
        best_nmse = np.inf
        for model in ['Logistic', 'Gompertz', 'Bertalanffy']:
            E = mouse_models[model]['nmse']
            i = E.index(min(E))
            best_params[mouse_id][model] = id2iterator[i]
            m = E[i]
            if m <= best_nmse:
                best_nmse = m
                best_model = model
        best_models[mouse_id] = best_model
    return best_params, best_models


def plot_tumor_data(all_data):
    volume_data = dict((m, all_data[m]['Volume']) for m in all_data.keys())
    df = pd.concat([pd.DataFrame(v, columns=[k]) for k, v in volume_data.items()], axis=1)
    df.plot(kind='line', legend=False, xlabel="Time point", ylabel="Tumor volume ")
    plt.show()


def save_all_model_fits(data):
    path = os.path.join(cwd, "model_params.json")
    with open(path, "w") as outfile: 
        json.dump(data, outfile)

def save_best_model_fits(all_models, model2iterator, mouse2model):
    iterators = {'fe': 0, 'lf': 1, 'be': 2}
    path =  os.path.join(cwd, "best_model_fits.txt")
    with open(path, "w") as f:
        f.write("mouse id \t model \t iterator \t params (a,b) \t nmse \n")
        for mouse_id in all_models.keys():
            id = mouse_id
            model = mouse2model[id]
            it = model2iterator[id][model]
            params = all_models[id][model]['params'][iterators[it]]
            e = all_models[id][model]['nmse'][iterators[it]]
            row = id + "\t" + model + "\t" + it + "\t" + str(params) + "\t" + str(e) + "\n"
            f.write(row)

def main():
    t = TumorGrowth()
    #plot_tumor_data(t.mouse2data)
    D1 = convert_all_models(t.mouse2model)
    D2, D3 = get_best_models(t.mouse2model)
    save_all_model_fits({'All Model Fits': D1, 'Best Iterators per Model': D2, 'Best Models': D3})
    save_best_model_fits(D1, D2, D3)
    for mouse in t.mouse2model.keys():
        t.plot(mouse, D2[mouse]['Logistic'], D2[mouse]['Gompertz'], D2[mouse]['Bertalanffy'])
    

if __name__ == "__main__":
    main()