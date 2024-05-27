import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import MaxNLocator

class sensitivity:
    def __init__(self, score_matrix, weights_arr, option_names, criteria_names):
        self.n_criteria = score_matrix.shape[0]
        self.n_options = score_matrix.shape[1]
        self.scores = score_matrix
        self.weights = weights_arr

        self.opt_names = option_names
        self.crit_names = criteria_names


    def weighted_score(self, option, weights):
        ws = np.sum(weights * option)
        return ws

    def perform_sensitivity_per_crit(self, criterion_pChanges, n_points: int = 50, plot_single=True, plot_summary=True):
        'initial weights'
        initial_weights = self.weights
        if n_points % 2 ==0:
            n_points +=1

        winner_difference_matrix = np.zeros((self.n_criteria, n_points))

        for i in range(self.n_criteria):
            indices = np.linspace(0, self.n_criteria-1, num=self.n_criteria, dtype=int)

            'initialise the final weights'
            final_weights = np.zeros((self.n_criteria, n_points))

            'get the ith % change and weight'
            percentage_change = criterion_pChanges[i]
            current_weight = initial_weights[i]

            'assign array of multiplied w_i to final weights'
            main_mult_range = np.linspace(1-percentage_change/100, 1+percentage_change/100, n_points)
            remaining_weights = initial_weights[indices != i]

            changed_main_w = main_mult_range * current_weight
            final_weights[i,:] = changed_main_w

            'assign array of multiplied w_k k!=i to final weights'
            remaining_mult_range = (1-changed_main_w)/(np.sum(remaining_weights))
            remaining_w = remaining_weights[:,np.newaxis] * remaining_mult_range
            final_weights[np.where(indices != i), :] = remaining_w

            'matmul over criteria'
            weighted_scores = np.einsum('ij, jk->ik', self.scores.T, final_weights)
            if plot_single:
                self.plot_sensitivity(x_axis=main_mult_range-1, weighted_scores=weighted_scores, idx=i, changes=criterion_pChanges)

            nominal_scores = weighted_scores[:, main_mult_range==1.].reshape(-1)
            nominal_difference = np.sort(nominal_scores)[-1] - np.sort(nominal_scores)[-2]

            nominal_winner_row = np.argwhere(nominal_scores == max(nominal_scores))[0,0]
            next_highest = np.max(weighted_scores[np.arange(self.n_options)!=nominal_winner_row, :], axis=0)
            differences = weighted_scores[nominal_winner_row, :] - next_highest

            p_differences = differences/weighted_scores[nominal_winner_row, :]*100
            winner_difference_matrix[i,:] = p_differences



        if plot_summary:
            y_labels = []
            for i, crit_name in enumerate(self.crit_names):
                y_labels.append(f'{crit_name}: ({self.weights[i]}'+r'$\pm$'+f'{criterion_pChanges[i]}%)')
            self.plot_summary_table(matrix=winner_difference_matrix, x_labels=main_mult_range-1, y_labels=y_labels)
        return initial_weights

    def plot_summary_table(self, matrix, x_labels = None, y_labels=None):
        print(matrix.shape, x_labels.shape)
        fig, ax = plt.subplots(1, 1, figsize=(10, 4), sharex=True, layout='constrained')

        #range = max(abs(np.min(matrix)), abs(np.max(matrix)))
        #divnorm=colors.TwoSlopeNorm(vcenter=0.,vmin=-range, vmax=range)
        abs_max = np.max(np.abs(matrix))
        rel_max = np.max([abs_max, 0])  # Positive relative max
        rel_min = -rel_max  # Negative relative min

        # Normalize the color map using the relative max and min values
        divnorm = colors.TwoSlopeNorm(vcenter=0.0, vmin=rel_min, vmax=rel_max)

        ims = ax.imshow(matrix, interpolation=None, cmap='PRGn', norm=divnorm, aspect='auto')

        cbar_ticks = np.array([np.min(matrix),0, np.max(matrix)])
        cbar = fig.colorbar(ims, ax=ax, orientation='horizontal', ticks = cbar_ticks)
        cbar.set_label('% change nominal difference best two options')


        ax.set_xlabel('Relative change')
        AR0 = 3
        if x_labels is not None:
            ax.set_xticks(np.arange(len(x_labels)))
            ax.xaxis.set_major_locator(MaxNLocator(integer=True, prune='both'))
            ax.set_xticklabels([f'{label:.1f}' for label in x_labels])

        if y_labels is not None:
            ax.set_yticks(np.arange(len(y_labels)))
            ax.set_yticklabels(y_labels)
        ax.set_aspect(1/ (AR0*matrix.shape[0]/matrix.shape[1]))
        plt.show()
        return

    def plot_sensitivity(self, idx, x_axis, weighted_scores, changes):
        #print(x_axis.shape, weighted_scores.shape, changes.shape)

        plt.figure(figsize=(5.5,5))
        for k in range(weighted_scores.shape[0]):
            plt.plot(x_axis, weighted_scores[k,:], label = f'{self.opt_names[k]}', linewidth=1)
        plt.title(f'Criterion: {self.crit_names[idx]}, change: {changes[idx]}%')
        plt.xlabel('relative change', fontsize=12)
        plt.ylabel('weighted score', fontsize=12)
        plt.legend(loc='lower left')
        plt.show()






def sens_structures():
    criterion_weights = np.array([.5, .2, .3])

    scores = np.array([[2, 3, 1,],
                       [4, 3, 2,],
                       [5, 1, 2,]])

    criterion_changes = np.array([100, 100, 100])

    design_option_names = ['truss+tower', 'truss+platform', 'branching', ]
    criteria_names = ['cost', 'maintenance', 'complexity']

    TEST = sensitivity(score_matrix=scores, weights_arr=criterion_weights, option_names=design_option_names,
                       criteria_names=criteria_names)
    TEST.perform_sensitivity_per_crit(criterion_pChanges=criterion_changes)


def sens_afc():
    criterion_weights = np.array([0.3, 0.3, 0.15, 0.15, 0.1])

    scores = np.array([[2, 2, 5, 5, 3 ,5, 5],
                    [4, 1, 1, 5 ,2, 1, 5],
                    [3, 4, 4, 3, 5, 2, 3 ],
                    [4, 3, 2, 4, 3, 1, 4],
                    [3, 4, 5, 5, 5, 2 ,2]])

    criterion_changes = np.array([100, 100, 100, 100, 100])
    design_option_names = ['single rotor', 'full sturcture', 'fixed HLD', 'retractable HLD', 'no HLD', 'fixed softwing', 'retractable softwing']
    criteria_names = ['vertical flow displacement', 'meteorologica versatility', 'reliability/complexity', 'structural resilience', 'innovation maturity']

    TEST = sensitivity(score_matrix=scores, weights_arr=criterion_weights, option_names=design_option_names, criteria_names=criteria_names)
    TEST.perform_sensitivity_per_crit(criterion_pChanges=criterion_changes)


def sens_rotor_types():
    criterion_weights = np.array([0.3, 0.1, 0.15, 0.15, 0.15, 0.15])

    scores = np.array([[2, 4, 3, 3],
                    [4, 4, 3, 3],
                    [3, 3, 4, 3],
                    [2, 2, 3 ,2],
                    [2, 2, 5, 3],
                    [2, 3 ,4, 2]])

    criterion_changes = np.array([50, 30, 70, 70, 60, 60])
    design_option_names = ['staggered', 'co-axial' , 'co-planar', 'wind wall']
    criteria_names = ['power generation', 'area efficiency', 'durability', 'mass', 'complexity', 'manufacturability']

    TEST = sensitivity(score_matrix=scores, weights_arr=criterion_weights, option_names=design_option_names, criteria_names=criteria_names)
    TEST.perform_sensitivity_per_crit(criterion_pChanges=criterion_changes)


def sens_rotor_number():
    criterion_weights = np.array([0.3, 0.1, 0.15, 0.15, 0.15, 0.15])
    #Original scores
    # scores = np.array([[3, 4, 5, 5],
    #                 [3, 4, 4, 5],
    #                 [2, 3, 4, 4],
    #                 [3, 4, 4, 3],
    #                 [4, 3, 3, 2],
    #                 [2, 3, 4, 4]])
   # updated scores 
    scores = np.array([[3, 4, 5, 5],
                    [3, 4, 4, 5],
                    [3, 4, 4, 5],
                    [3, 4, 4, 3],
                    [4, 4, 3, 3],
                    [2, 3, 4, 5]])

    criterion_changes = np.array([30, 30, 50, 60, 50, 40])
    design_option_names = ['11 rotors', '23 rotors' , '33 rotors', '53 rotors']
    criteria_names = ['power generation', 'area efficiency', 'durability', 'mass', 'complexity', 'manufacturability']

    TEST = sensitivity(score_matrix=scores, weights_arr=criterion_weights, option_names=design_option_names, criteria_names=criteria_names)
    TEST.perform_sensitivity_per_crit(criterion_pChanges=criterion_changes)


def sens_yaw_control():
    criterion_weights = np.array([0.2, 0.3, 0.15, 0.2, 0.15])

    scores = np.array([[1, 4, 3, 2, 4, 3],
                       [2, 3, 4, 5, 3, 4],
                       [3, 4, 3, 4, 4, 3],
                       [5, 3, 4, 5, 3, 4],
                       [2, 3, 1, 3, 4, 1]])

    # good:
    # bad: failure rate, power req, versatility

    criterion_changes = np.array([100, 50, 100, 100, 50])
    design_option_names = ['bearing motor', 'bearing differential pitch', 'bearing reverse thrust', 'turntable motor', 'turntable differential pitch', 'turntable reverse thrust']
    criteria_names = ['power required', 'response time', 'failure rate', 'versatility', 'complexity']

    TEST = sensitivity(score_matrix=scores, weights_arr=criterion_weights, option_names=design_option_names, criteria_names=criteria_names)
    TEST.perform_sensitivity_per_crit(criterion_pChanges=criterion_changes)


def sense_generator():
    criterion_weights = np.array([0.1, 0.1, 0.15, 0.2, 0.15, 0.15, 0.1, 0.05])

    scores = np.array([[2, 2, 5, 5, 1],
                       [5, 5, 4, 4, 1],
                       [3, 3, 4, 4, 2],
                       [4, 3, 4, 4, 5],
                       [5, 4, 3, 3, 2],
                       [3, 4, 2, 1, 5],
                       [4, 4, 4, 1, 4],
                       [4, 4, 3, 1, 3]])
    
    criterion_changes = np.array([50, 50, 100, 100, 50, 50, 50, 50])
    design_option_names = ['brushless DFIG', 'DFIG', 'squirrel cage induction generator', 'permanent magnet synchronous generator', 'constant speed SCIG']
    criteria_names = ['speed variation', 'power factor and voltage conrol', 'energy extraction', 'maintainability and reliability', 'efficiency losses', 'cost', 'mass', 'sustainability']
    # good: speed var, PF,  , eff, costs, mass, sust
    # BAD: energy extraction, maint,

    TEST = sensitivity(score_matrix=scores, weights_arr=criterion_weights, option_names=design_option_names, criteria_names=criteria_names)
    TEST.perform_sensitivity_per_crit(criterion_pChanges=criterion_changes)


def sense_pitch():
    criterion_weights = np.array([0.3, 0.35, 0.2, 0.15])

    scores = np.array([[5, 4],
                       [2, 4],
                       [2, 4],
                       [1, 3]])

    criterion_changes = np.array([50, 50, 50, 50])
    design_option_names = ['individual pitch control', 'group pitch control']
    criteria_names = ['blade fatigue', 'cost', 'reliability', 'complexity']

    TEST = sensitivity(score_matrix=scores, weights_arr=criterion_weights, option_names=design_option_names, criteria_names=criteria_names)
    TEST.perform_sensitivity_per_crit(criterion_pChanges=criterion_changes)


def sense_drive_train():
    criterion_weights = np.array([0.1, 0.15, 0.1, 0.2, 0.15, 0.1, 0.1, 0.1])

    scores = np.array([[4, 3, 5, 5],
                       [1, 1, 4, 4],
                       [2, 3, 4, 4],
                       [5, 3, 4, 4],
                       [4, 4, 3, 3],
                       [4, 3, 2, 1],
                       [2, 3, 4, 1],
                       [4, 3, 3, 1]])

    criterion_changes = np.array([50, 100, 50, 50, 50, 50, 100, 100])
    design_option_names = ['modular', 'partially integrated', 'integrated', 'tbd']
    criteria_names = ['torque isolation', 'complexity', 'alignment tolerances', 'maintainance', 'interdependencies', 'cost', 'mass', 'sustainability']
    # BAD: complex, mass, sustainability,

    # GOOD: ti, alingment, inter, maint, cost
    TEST = sensitivity(score_matrix=scores, weights_arr=criterion_weights, option_names=design_option_names, criteria_names=criteria_names)
    TEST.perform_sensitivity_per_crit(criterion_pChanges=criterion_changes)


def system_trade_off(score_change = False):
    criterion_weights = np.array([0.3, 0.25, 0.2, 0.25])

    if not score_change:
        scores = np.array([[4,5,2],
                           [3,4,1],
                           [3,4,2],
                           [2,2,4]])
    else:
        scores = np.array([[3, 4, 2],
                           [3, 4, 1],
                           [3, 3, 2],
                           [2, 2, 3]])


    criterion_changes = np.array([100, 50, 50, 50])
    design_option_names = ['in-plane', 'out-plane', 'platform/grid/turntable/motor']
    criteria_names = ['cost', 'complexity', 'sustainability', 'maintainability',]

    TEST = sensitivity(score_matrix=scores, weights_arr=criterion_weights, option_names=design_option_names, criteria_names=criteria_names)
    TEST.perform_sensitivity_per_crit(criterion_pChanges=criterion_changes)



if __name__ == '__main__':
    #sens_structures()
    #system_trade_off(score_change=False)
    sens_rotor_types()
    #sense_drive_train()
    #sense_generator()

    #sense_pitch()

    #sense_pitch()

    sens_rotor_types()