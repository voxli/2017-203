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

    double temp1, temp2;

    for (int i = 0; i < N; i++)
    {
        temp1 = 0;
        temp2 = 1;

        for (int k = 0; k < i; k++)
        {
            temp1 += pow(A[k][i], 2);
        }

        A[i][i] = sqrt(A[i][i] - double(temp1));

        if(i == 0)
        {
            temp2 = 0;
        }
        else
        {
            for(int l = i; l < N; l++)
            {
                temp2 = temp2 * A[i - 1][l];
            }
        }

        for (int j = 0; j < N; j++)
        {
            if (j < i)
            {
                A[i][j] = 0;
            }
            else if (i == j)
            {
                continue;
            }
            else
            {
                A[i][j] = (A[i][j] - temp2) / A[i][i];
            }
        }
    }

    for (int i = 0; i < N; i++)
    {
        temp1 = 0;

        for (int k = 0; k < i; k++)
        {
            temp1 = temp1 + A[k][i] * b[k];
        }

        b[i] = (b[i] - temp1) / A[i][i];
    }

    for (int k = N - 1; k >= 0; k--)
    {
        double resault = 0;

        for (int i = k + 1; i < N; i++)
        {
            resault += A[k][i] * x[i];
        }

        x[k] = (b[k] - resault) / A[k][k];
    }

}



/**
 * Метод прогонки
 */
void kidyaevvs::lab4()
{

    double* new_A = new double[N];
    double* new_b = new double[N];

    new_A[0] = A[0][1] / (-A[0][0]);
    new_b[0] = b[0] / A[0][0];

    for(int i = 1; i < N; i++)
    {
        new_A[i] = A[i][i+1] / (-A[i][i-1] * new_A[i-1] - A[i][i]);
        new_b[i] = (-b[i] + A[i][i-1] * new_b[i-1]) / ( -A[i][i-1] * new_A[i-1] - A[i][i]);
    }

    for(int i = N - 1; i >= 0; i--)
    {
        x[i] = new_A[i] * x[i+1] + new_b[i];
    }

    delete[] new_A;
    delete[] new_b;

}


/**
 * Метод Якоби
 */
void kidyaevvs::lab5()
{
    double Eps = 0.0001;
    double* temp_x = new double[N];
    double norm;

    for (int i = 0; i < N; i++)
    {
        x[i]=0;
    }

    do {
        for (int i = 0; i < N; i++)
        {
            temp_x[i] = b[i];

            for (int j = 0; j < N; j++)
            {
                if (i != j)
                {
                    temp_x[i] -= A[i][j] * x[j];
                }
            }

            temp_x[i] /= A[i][i];
        }

        norm = fabs(x[0] - temp_x[0]);

        for (int i = 0; i < N; i++)
        {
            if (fabs(x[i] - temp_x[i]) > norm)
            {
                norm = fabs(x[i] - temp_x[i]);
            }

            x[i] = temp_x[i];
        }
    } while (norm > Eps);

    delete[] temp_x;
}



/**
 * Метод Зейделя
 */
void kidyaevvs::lab6()
{
    double eps = 0.0001;
    double* y = new double[N];
    double norm = 0;
    double temp = 0;

    for (int i = 0; i < N; i++)
    {
        x[i] = 0;
    }

    do
    {
        for (int i = 0; i < N; i++)
        {
            y[i] = x[i];
        }

        for (int i = 0; i < N; i++)
        {
            temp = 0;
            norm = 0;

            for (int j = 0; j < i; j++)
            {
                temp += (A[i][j] * x[j]);
            }

            for (int j = i + 1; j < N; j++)
            {
                temp += (A[i][j] * x[j]);
            }

            x[i] = (b[i] - temp) / A[i][i];

            for (int i = 0; i < N; i++)
            {
                norm += (x[i] - y[i])*(x[i] - y[i]);
            }
        }
    } while (sqrt(norm) >= eps);

    delete[] y;
}



/**
 * Один из градиентных методов
 */
void kidyaevvs::lab7()
{

}

/**
 * Один из градиентных методов
 */
void kidyaevvs::lab8()
{

}
