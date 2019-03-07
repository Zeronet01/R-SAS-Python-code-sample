import SimPy.RandomVariantGenerators as RVGs
import SimPy.MarkovClasses as Markov
import SimPy.SamplePathClasses as PathCls
from InputData import HealthState


class Patient:
    def __init__(self, id, trans_rate_matrix):
        """ initiates a patient
        :param id: ID of the patient
        :param trans_rate_matrix: transition rate matrix
        """
        self.id = id
        self.rng = RVGs.RNG(seed=id)  # random number generator for this patient
        # gillespie algorithm
        self.gillespie = Markov.Gillespie(transition_rate_matrix=trans_rate_matrix)
        self.stateMonitor = PatientStateMonitor()  # patient state monitor

    def simulate(self, sim_length):
        """ simulate the patient over the specified simulation length """

        t = 0  # simulation time
        if_stop = False  # set to true to stop the simulation

        # while the patient is alive and simulation length is not yet reached
        while not if_stop:

            dt, new_state_index = self.gillespie.get_next_state(
                current_state_index=self.stateMonitor.currentState.value,
                rng=self.rng)

            if dt is None or dt + t > sim_length:
                if_stop = True
            else:
                # increment time
                t += dt
                # update health state
                self.stateMonitor.update(time=t, new_state=HealthState(new_state_index))


class PatientStateMonitor:
    """ to update patient outcomes (years survived, cost, etc.) throughout the simulation """
    def __init__(self):

        self.currentState = HealthState.Well    # current health state
        self.survivalTime = None      # survival time
        self.timeToSTROKES = None        # time to develop STROKES
        self.ifDevelopedSTROKES = False  # if the patient developed STROKES
        self.num_strokes = 0

    def update(self, time, new_state):
        """
        update the current health state to the new health state
        :param time: current time
        :param new_state: new state
        """

        # update survival time
        if new_state == HealthState.Stroke_Death or HealthState.NATUAL_DEATH:
            self.survivalTime = time

        # update time until STROKES
        if self.currentState != HealthState.Stroke and new_state == HealthState.Stroke:
            self.ifDevelopedSTROKES = True
            self.timeToSTROKES = time

        # update current health state
        self.currentState = new_state

        # update times of strokes
        if new_state == HealthState.Stroke:
            self.num_strokes +=1


class Cohort:
    def __init__(self, id, pop_size, transition_matrix):
        """ create a cohort of patients
        :param id: cohort ID
        :param pop_size: population size of this cohort
        :param transition_matrix: probability transition matrix
        """
        self.id = id
        self.patients = []  # list of patients
        self.cohortOutcomes = CohortOutcomes()  # outcomes of the this simulated cohort

        # populate the cohort
        for i in range(pop_size):
            # create a new patient (use id * pop_size + n as patient id)
            patient = Patient(id=id * pop_size + i, trans_rate_matrix=transition_matrix)
            # add the patient to the cohort
            self.patients.append(patient)

    def simulate(self, sim_length):
        """ simulate the cohort of patients over the specified number of time-steps
        :param sim_length: simulation length
        """
        # simulate all patients
        for patient in self.patients:
            # simulate
            patient.simulate(sim_length)

        # store outputs of this simulation
        self.cohortOutcomes.extract_outcomes(self.patients)


class CohortOutcomes:
    def __init__(self):

        self.survivalTimes = []         # patients' survival times
        self.timesToSTROKES = []           # patients' times to STROKES
        self.meanSurvivalTime = None    # mean survival times
        self.meanTimeToSTROKES = None      # mean time to STROKES
        self.nLivingPatients = None     # survival curve (sample path of number of alive patients over time)
        self.Num_stroke=[]          # patients' number of stokes
        self.MeanNumberSTROKES = None      # mean number of STROKES

    def extract_outcomes(self, simulated_patients):
        """ extracts outcomes of a simulated cohort
        :param simulated_patients: a list of simulated patients"""

        # record survival time and time until STROKES
        for patient in simulated_patients:
            if not (patient.stateMonitor.survivalTime is None):
                self.survivalTimes.append(patient.stateMonitor.survivalTime)
            if patient.stateMonitor.ifDevelopedSTROKES:
                self.timesToSTROKES.append(patient.stateMonitor.timeToSTROKES)
            if patient.stateMonitor.ifDevelopedSTROKES:
                self.Num_stroke.append(patient.stateMonitor.num_strokes)

        # calculate mean survival time
        self.meanSurvivalTime = sum(self.survivalTimes) / len(self.survivalTimes)
        # calculate mean time to STROKES
        self.meanTimeToSTROKES = sum(self.timesToSTROKES)/len(self.timesToSTROKES)
        # calculate mean number of STROKES
        self.MeanNumberSTROKES=sum(self.Num_stroke)/len(self.Num_stroke)

        # survival curve
        self.nLivingPatients = PathCls.PrevalencePathBatchUpdate(
            name='# of living patients',
            initial_size=len(simulated_patients),
            times_of_changes=self.survivalTimes,
            increments=[-1]*len(self.survivalTimes)
        )
