const int N_SpeechGoals = 8;
const char SpeechGoalName[N_SpeechGoals] = {'i', 'e', '3', 'a', '@', 'O', 'k', 't'}; 

void ellipticzone_error(char target_code, char target_modality, double *y, double *err, double scaling=1.0);
void ellipticzone_error_ds(char target_code, char target_modality, double *y, double **H, double **He, double scaling=1.0);
void init_goal_geometry();