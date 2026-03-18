#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif


//измерение времени
double get_precise_time_ms() {
#ifdef _WIN32
    LARGE_INTEGER frequency, counter;
    QueryPerformanceFrequency(&frequency);
    QueryPerformanceCounter(&counter);
    return (double)counter.QuadPart * 1000.0 / (double)frequency.QuadPart;
#else
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec * 1000.0 + (double)tv.tv_usec / 1000.0;
#endif
}

//вычисление нормы вектора(для матрицы Гильберта)
double vector_norm(double* x, int n) {
    double sum = 0;
    for (int i = 0; i < n; i++) {
        sum += x[i] * x[i];
    }
    return sqrt(sum);
}

//создание матрицы в динамической памяти
double** create_matrix(int n) {
    double** mat = (double**)malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) {
        mat[i] = (double*)malloc(n * sizeof(double));
    }
    return mat;
}

//создание вектора
double* create_vector(int n) {
    return (double*)malloc(n * sizeof(double));
}

//освобождение памяти матрицы
void free_matrix(double** mat, int n) {
    for (int i = 0; i < n; i++) {
        free(mat[i]);
    }
    free(mat);
}

//освобождение вектора
void free_vector(double* vec) {
    free(vec);
}

//копирование матрицы
void copy_matrix(double** src, double** dst, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            dst[i][j] = src[i][j];
        }
    }
}

//копирование вектора
void copy_vector(double* src, double* dst, int n) {
    for (int i = 0; i < n; i++) {
        dst[i] = src[i];
    }
}

//генерация случайной матрицы с элементами из [-1, 1]
void generate_random_matrix(double** A, int n, unsigned int seed) {
    srand(seed);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;  
        }
    }
}

//генерация случайного вектора
void generate_random_vector(double* b, int n, unsigned int seed) {
    srand(seed + 1000);  
    for (int i = 0; i < n; i++) {
        b[i] = ((double)rand() / RAND_MAX) * 2.0 - 1.0; 
    }
}

//генерация матрицы Гильберта 
void generate_hilbert_matrix(double** H, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            H[i][j] = 1.0 / (i + j + 1);
        }
    }
}

//умножение матрицы на вектор (для проверки невязки)
void matrix_vector_mult(double** A, double* x, double* res, int n) {
    for (int i = 0; i < n; i++) {
        res[i] = 0;
        for (int j = 0; j < n; j++) {
            res[i] += A[i][j] * x[j];
        }
    }
}

//вычисление невязки
double compute_residual(double** A, double* x, double* b, int n) {
    double* Ax = create_vector(n);
    matrix_vector_mult(A, x, Ax, n);
    
    double max_res = 0;
    for (int i = 0; i < n; i++) {
        double res = fabs(Ax[i] - b[i]);
        if (res > max_res) max_res = res;
    }
    
    free_vector(Ax);
    return max_res;
}

//LU-разложение (без перестановок)
int lu_decomposition(double** A, double** L, double** U, int n) {
    //инициализация
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            L[i][j] = 0;
            U[i][j] = 0;
        }
        L[i][i] = 1;  
    }
    
    //основной алгоритм
    for (int k = 0; k < n; k++) {
        //вычисление U для строки k
        for (int j = k; j < n; j++) {
            double sum = 0;
            for (int s = 0; s < k; s++) {
                sum += L[k][s] * U[s][j];
            }
            U[k][j] = A[k][j] - sum;
        }
        
        //проверка на ноль
        if (fabs(U[k][k]) < 1e-15) {
            printf("Ошибка: матрица требует перестановок\n");
            return 0;
        }
        
        //вычисление L для столбца k
        for (int i = k + 1; i < n; i++) {
            double sum = 0;
            for (int s = 0; s < k; s++) {
                sum += L[i][s] * U[s][k];
            }
            L[i][k] = (A[i][k] - sum) / U[k][k];
        }
    }
    
    return 1;
}

//прямая подстановка (Ly = b)
void forward_substitution(double** L, double* b, double* y, int n) {
    for (int i = 0; i < n; i++) {
        double sum = 0;
        for (int j = 0; j < i; j++) {
            sum += L[i][j] * y[j];
        }
        y[i] = b[i] - sum;  //L[i][i] = 1
    }
}

//обратная подстановка (Ux = y)
void backward_substitution(double** U, double* y, double* x, int n) {
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0;
        for (int j = i + 1; j < n; j++) {
            sum += U[i][j] * x[j];
        }
        x[i] = (y[i] - sum) / U[i][i];
    }
}

//решение системы через LU-разложение
void solve_lu(double** L, double** U, double* b, double* x, int n) {
    double* y = create_vector(n);
    forward_substitution(L, b, y, n);
    backward_substitution(U, y, x, n);
    free_vector(y);
}

//печать таблицы результатов
void print_results_table(const char* title, int* sizes, int num_sizes,
                         double** lu_decomp_times, double** lu_solve_times,
                         double** lu_total_times) {
    printf("\n%s\n", title);
    printf("|  n  | LU-разложение | LU-решение | LU-сумма |\n");
    printf("--------------------------------------------------\n");
    
    for (int i = 0; i < num_sizes; i++) {
        printf("| %3d | %13.3f | %10.3f | %9.3f |\n",
               sizes[i], lu_decomp_times[i][0], lu_solve_times[i][0], lu_total_times[i][0]);
    }
    
}

//экс-т 1: сравнение времени для одной системы (только LU)
void experiment_1() {
    printf("\n Эксперимент 1: Время LU-разложения для одной системы \n");
    
    int sizes[] = {100, 200, 500};  
    int num_sizes = sizeof(sizes) / sizeof(sizes[0]);
    unsigned int seed = 12345; 
    
    //массивы для хранения результатов
    double* lu_decomp_times[3];
    double* lu_solve_times[3];
    double* lu_total_times[3];
    
    for (int i = 0; i < num_sizes; i++) {
        lu_decomp_times[i] = create_vector(1);
        lu_solve_times[i] = create_vector(1);
        lu_total_times[i] = create_vector(1);
    }
    
    for (int idx = 0; idx < num_sizes; idx++) {
        int n = sizes[idx];
        printf("\nТестирование n = %d...\n", n);
        
        double** A = create_matrix(n);
        double** L = create_matrix(n);
        double** U = create_matrix(n);
        double* b = create_vector(n);
        double* x = create_vector(n);
        
        generate_random_matrix(A, n, seed);
        generate_random_vector(b, n, seed);
        
        //LU-разложение + решение
        double start = get_precise_time_ms();
        
        double decomp_start = get_precise_time_ms();
        if (!lu_decomposition(A, L, U, n)) {
            printf("  Ошибка LU-разложения\n");
            lu_decomp_times[idx][0] = -1;
            lu_solve_times[idx][0] = -1;
            lu_total_times[idx][0] = -1;
            
            free_matrix(A, n);
            free_matrix(L, n);
            free_matrix(U, n);
            free_vector(b);
            free_vector(x);
            continue;
        }
        double decomp_end = get_precise_time_ms();
        lu_decomp_times[idx][0] = decomp_end - decomp_start;
        
        double solve_start = get_precise_time_ms();
        solve_lu(L, U, b, x, n);
        double solve_end = get_precise_time_ms();
        lu_solve_times[idx][0] = solve_end - solve_start;
        
        double end = get_precise_time_ms();
        lu_total_times[idx][0] = end - start;
        
        printf("  LU-разложение: %.3f мс\n", lu_decomp_times[idx][0]);
        printf("  LU-решение: %.3f мс\n", lu_solve_times[idx][0]);
        printf("  LU-сумма: %.3f мс\n", lu_total_times[idx][0]);
        
        //проверка невязки
        double resid = compute_residual(A, x, b, n);
        printf("  Невязка: %.2e\n", resid);
        
        free_matrix(A, n);
        free_matrix(L, n);
        free_matrix(U, n);
        free_vector(b);
        free_vector(x);
    }
    
    //вывод таблицы
    print_results_table("РЕЗУЛЬТАТЫ Эксперимента 1", sizes, num_sizes,
                        lu_decomp_times, lu_solve_times, lu_total_times);
    
    for (int i = 0; i < num_sizes; i++) {
        free_vector(lu_decomp_times[i]);
        free_vector(lu_solve_times[i]);
        free_vector(lu_total_times[i]);
    }
}

//экс-т 2: множественные правые части (только LU)
void experiment_2() {
    printf("\n Эксперимент 2: Множественные правые части (LU-разложение) \n");
    
    int n = 500;  
    int k_vals[] = {1, 10, 100};
    int num_k = sizeof(k_vals) / sizeof(k_vals[0]);
    unsigned int seed = 54321;
    
    printf("\nРазмер матрицы: %d\n", n);
    printf("|  k  | LU-разложение | LU-решение | Сумма |\n");
    printf("---------------------------------------------\n");
     
    double** A = create_matrix(n);
    generate_random_matrix(A, n, seed);
    
    //LU-разложение делаем один раз
    double** L = create_matrix(n);
    double** U = create_matrix(n);
    double lu_decomp_time;
    
    double start = get_precise_time_ms();
    if (!lu_decomposition(A, L, U, n)) {
        printf("Ошибка LU-разложения\n");
        free_matrix(A, n);
        free_matrix(L, n);
        free_matrix(U, n);
        return;
    }
    double end = get_precise_time_ms();
    lu_decomp_time = end - start;
    printf("Время однократного LU-разложения: %.3f мс\n\n", lu_decomp_time);
    
    for (int k_idx = 0; k_idx < num_k; k_idx++) {
        int k = k_vals[k_idx];
        
        //генерация k правых частей
        double** b_list = create_matrix(k);  // k x n
        double** x_results = create_matrix(k);  
        
        for (int i = 0; i < k; i++) {
            b_list[i] = create_vector(n);
            x_results[i] = create_vector(n);
            generate_random_vector(b_list[i], n, seed + i);
        }
        
        //Решение через LU (только подстановки для каждой правой части)
        double lu_solve_start = get_precise_time_ms();
        for (int i = 0; i < k; i++) {
            solve_lu(L, U, b_list[i], x_results[i], n);
        }
        double lu_solve_end = get_precise_time_ms();
        double lu_solve_total = lu_solve_end - lu_solve_start;
        double lu_total = lu_decomp_time + lu_solve_total;
        
        printf("| %3d | %13.3f | %11.3f | %7.3f |\n",
               k, lu_decomp_time, lu_solve_total, lu_total);
        
        for (int i = 0; i < k; i++) {
            free_vector(b_list[i]);
            free_vector(x_results[i]);
        }
        free_matrix(b_list, k);
        free_matrix(x_results, k);
    }
    
    printf("---------------------------------------------\n");
    
    free_matrix(A, n);
    free_matrix(L, n);
    free_matrix(U, n);
}

//экс-т 3: проверка точности на матрице Гильберта (только LU)
void experiment_3() {
    printf("\n Эксперимент 3: Проверка точности на матрице Гильберта \n");
    
    int sizes[] = {5, 10, 15};
    int num_sizes = sizeof(sizes) / sizeof(sizes[0]);

    printf("|  n  | Погрешность | Невязка     |\n");
    printf("--------------------------------------------\n");
    
    for (int idx = 0; idx < num_sizes; idx++) {
        int n = sizes[idx];
        
        //создание матрицы Гильберта
        double** H = create_matrix(n);
        generate_hilbert_matrix(H, n);
        
        double* x_exact = create_vector(n);
        for (int i = 0; i < n; i++) {
            x_exact[i] = 1.0;
        }
        
        //вычисление правой части b = H * x_exact
        double* b = create_vector(n);
        matrix_vector_mult(H, x_exact, b, n);
        
        double* x = create_vector(n);
        
        //LU-разложение
        double** L = create_matrix(n);
        double** U = create_matrix(n);
        
        if (lu_decomposition(H, L, U, n)) {
            solve_lu(L, U, b, x, n);
            
            //вычисление погрешности
            double* diff = create_vector(n);
            for (int i = 0; i < n; i++) {
                diff[i] = x[i] - x_exact[i];
            }
            double error = vector_norm(diff, n) / vector_norm(x_exact, n);
            
            //вычисление невязки
            double resid = compute_residual(H, x, b, n);
            
            printf("| %3d | %11.2e | %11.2e |\n", n, error, resid);
            
            free_vector(diff);
        } else {
            printf("| %3d |   ошибка    |   ошибка    |\n", n);
        }
        
        printf("--------------------------------------------\n");
        
        free_matrix(H, n);
        free_matrix(L, n);
        free_matrix(U, n);
        free_vector(x_exact);
        free_vector(b);
        free_vector(x);
    }
}

// Основная функция
int main() {
    printf("  ЛАБОРАТОРНАЯ РАБОТА: LU-РАЗЛОЖЕНИЕ\n");
    
    experiment_1();
    
    experiment_2();
    
    experiment_3();
    
    printf("\nРабота завершена.\n");
    
    
    return 0;
}
