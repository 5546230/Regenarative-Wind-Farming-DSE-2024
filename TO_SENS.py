import numpy as np
import matplotlib.pyplot as plt


criterion_weights = np.array([.5, .2, .3])

scores = np.array([[3, 2, 1, 5],
                   [4, 3, 2, 1],
                   [5, 1, 2, 4]])

criterion_changes = np.array([50, 50, 50])
design_option_names = ['a', 'b', 'c', 'd']
criteria_names = ['i', 'j', 'k']
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

    def perform_sensitivity(self, criterion_pChanges, n_points: int = 10):
        initial_weights = self.weights.copy()

        for i in range(self.n_criteria):
            indices = np.linspace(0, self.n_criteria-1, self.n_criteria, dtype=int)
            final_weights = np.zeros((self.n_criteria, n_points))

            percentage_change = criterion_pChanges[i]
            current_weight = initial_weights[i]


            first_mult_range = np.linspace(1-percentage_change/100, 1+percentage_change/100, n_points)
            remaining_weights = initial_weights[initial_weights!=current_weight]
            changed_w = first_mult_range * current_weight
            final_weights[i,:] = changed_w


            remaining_mult_range = (1-changed_w)/(np.sum(remaining_weights))
            remaining_w = remaining_weights[:,np.newaxis] * remaining_mult_range
            final_weights[indices[self.weights != current_weight],:] = remaining_w

            weighted_scores = np.einsum('ij, jk->ik', self.scores.T, final_weights)
            self.plot_sensitivity(x_axis=first_mult_range, weighted_scores=weighted_scores, idx=i)
        return

    def plot_sensitivity(self, idx, x_axis, weighted_scores):
        for k in range(weighted_scores.shape[0]):
            plt.plot(x_axis, weighted_scores[k,:], label = f'{self.opt_names[k]}')
        plt.title(f'{self.crit_names[idx]}')
        plt.xlabel('relative change')
        plt.ylabel('weighted score')
        plt.legend()
        plt.show()





TEST = sensitivity(score_matrix=scores, weights_arr=criterion_weights, option_names=design_option_names, criteria_names=criteria_names)
TEST.perform_sensitivity(criterion_pChanges=criterion_changes)
