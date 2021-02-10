// fmin.cpp
/*
			
		Файл содержит функцию fmin, задачей которой является поиск точки минимума
		функции одной переменной.Автор: Ненашев Олег. Группа 2081/2
*/

  #include<math.h>
   double fmin(double beg,					    // Начало интервала поиска 
			double end,							// Конец интервала поиска
			double (*f) (double x),				// Заданная функция
			double Tol,							// Погрешность
			double& Count,						// Значение функции в точке минимума
			int& flag)							// Флаг
{
	double	BegBuf=f(beg),						// Значение функции в начале интервала
			EndBuf=f(end),						// Значение функции в начале интервала
			Lcenter,Rcenter, Lfcenter,Rfcenter;	// Координаты точек
	int counter=0;								// Счётчик
	int maxcounts=500;							// Максимальное число итераций

	flag=0;




	// Поиск минимального значения на участке
	do
	{
		// Анализ середины интервала
		Lcenter=(end + beg - Tol)/2;
		Lfcenter=f(Lcenter);

		Rcenter=(end + beg + Tol)/2;
		Rfcenter=f(Rcenter);
		if(Lfcenter<Rfcenter)
		{
			end=Rcenter;
			EndBuf=Rfcenter;
		}

		// Выбор нового интервала
		else if(Lfcenter>Rfcenter)
		{
			beg=Lcenter;
			BegBuf=Lfcenter;
		}
		else
		{
			return (Lcenter+Rcenter)/2;
		}

		// Анализ условия выхода
		if(fabs(end-beg)<2*Tol)
		{
			Lcenter=(end+beg)/2;
			Count=f(Lcenter);
			return Lcenter;
		}	

		counter++;
	}while(counter<maxcounts);

	flag=1;
	return end;

}