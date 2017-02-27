#include "kidyaevvs.h"

/**
 * Метод Гаусса
 */
void kidyaevvs::lab1()
{
	double new_var;

	for (int k = 0; k < N; k++)
	{
		for(int i = k + 1; i < N; i++)
		{
			new_var = A[i][k] / A[k][k];
			b[i] -= b[k] * new_var;

			for(int j = 0; j < N; j++)
			{
				A[i][j] -= A[k][j] * new_var;
			}
		}
	}

	for(int i = N - 1; i >= 0; i--)
	{
		x[i] = b[i] / A[i][i];

		for(int j = i + 1; j < N; j++)
		{
			x[i] += -A[i][j] * x[j] / A[i][i];
		}
	}

}


/**
 * Метод Гаусса с выбором главного элемента
 */
void kidyaevvs::lab2()
{
	for (int k = 0; k < N; k++)
	{
		int id_str_of_max_el = -1;
        double max_el = 0;

        for (int p = 0; p < N; p++)
		{
            if (abs(A[k][p]) >= max_el)
			{
                max_el = abs(A[k][p]);
                id_str_of_max_el = p;
            }
        }

        if (id_str_of_max_el != -1)
		{
            for (int j = 0; j < N; j++)
			{
               std :: swap(A[j][id_str_of_max_el], A[j][k]);
            }

            std :: swap(b[id_str_of_max_el], b[k]);

            for (int i = k + 1; i < N; i++)
            {
                double coef = A[i][k] / (A[k][k] * 1.0);

                for (int j = k; j < N; j++)
                {
					A[i][j] -= coef * A[k][j];
				}
                b[i] -= coef * b[k];
            }

        }

    }


    x[N - 1] = b[N - 1];

    for (int i = N - 2; i >= 0; i--)
	{
        x[i] = b[i];

        for (int j = i + 1; j < N; j++)
        {
			x[i] -= x[j] * A[i][j];
		}
        x[i] /= A[i][i];
    }
}



/**
 * Метод квадратного корня (метод Холецкого)
 */
void kidyaevvs::lab3()
{

}



/**
 * Метод прогонки
 */
void kidyaevvs::lab4()
{

}



/**
 * Метод Якоби
 */
void kidyaevvs::lab5()
{

}



/**
 * Метод Зейделя
 */
void kidyaevvs::lab6()
{

}



/**
 * Один из градиентных методов
 */
void kidyaevvs::lab7()
{

}
