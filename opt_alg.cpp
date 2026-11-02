#include"opt_alg.h"

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	// Zmienne wej�ciowe:
	// ff - wska�nik do funkcji celu
	// N - liczba zmiennych funkcji celu
	// lb, ub - dolne i g�rne ograniczenie
	// epslion - zak��dana dok�adno�� rozwi�zania
	// Nmax - maksymalna liczba wywo�a� funkcji celu
	// ud1, ud2 - user data
	try
	{
		solution Xopt;
		Xopt = rand_mat(N);									// losujemy macierz Nx1 stosuj�c rozk�ad jednostajny na przedziale [0,1]
		for (int i = 0; i < N; ++i)							// przeskalowywujemy rozwi�zanie do przedzia�u [lb, ub]
			Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
		Xopt.fit_fun(ff, ud1, ud2);							// obliczamy warto�� funkcji celu
		Xopt.ud = (Xopt.y);									// warto�� funckji celu zapami�tujemy w macierzy ud
		if (Xopt.y < epsilon)								// sprawdzamy 1. kryterium stopu
		{
			Xopt.flag = 1;
			return Xopt;
		}
		solution X;
		while (true)
		{
			X = rand_mat(N);									// losujemy macierz Nx1 stosuj�c rozk�ad jednostajny na przedziale [0,1]
			for (int i = 0; i < N; ++i)							// przeskalowywujemy rozwi�zanie do przedzia�u [lb, ub]
				X.x(i) = (ub(i) - lb(i)) * X.x(i) + lb(i);
			X.fit_fun(ff, ud1, ud2);							// obliczmy warto�� funkcji celu
			if (X.y < Xopt.y)
			{
				Xopt = X;
				if (Xopt.y < epsilon)							// sprawdzmy 1. kryterium stopu
				{
					Xopt.flag = 1;								// flaga = 1 ozancza znalezienie rozwi�zanie z zadan� dok�adno�ci�
					break;
				}
			}
			if (solution::f_calls > Nmax)						// sprawdzmy 2. kryterium stopu
			{
				Xopt.flag = 0;									// flaga = 0 ozancza przekroczenie maksymalnej liczby wywo�a� funkcji celu
				break;
			}
			Xopt.ud.add_row(Xopt.y);							// warto�� funckji dodajemy jako kolejny wiersz macierzy ud
		}
		Xopt.ud.add_row(Xopt.y);
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}

double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		double* p = new double[2] { 0, 0 };
		solution X0, X1;
		X0.x = x0;
		X0.fit_fun(ff, ud1, ud2);
		X1.x = x0 + d;
		X1.fit_fun(ff, ud1, ud2);
		
		if (m2d(X1.y) == m2d(X0.y))
		{
			p[0] = m2d(X0.x);
			p[1] = m2d(X1.x);
			return p;
		}
		if (m2d(X1.y) > m2d(X0.y))
		{
			d = -d;
			X1.x = m2d(X0.x) + d;
			X1.fit_fun(ff, ud1, ud2);
			if (m2d(X1.y) >= m2d(X0.y))
			{
				p[0] = m2d(X1.x);
				p[1] = m2d(X0.x) - d;
				return p;
			}
		}
		
		int i = 1;
		solution X_prev = X0, X_curr = X1, X_next;
		while (true)
		{
			if (solution::f_calls > Nmax)
				throw string("Przekroczono Nmax");
				
			double x_val = x0 + pow(alpha, i) * d;
			X_next.x = x_val;
			X_next.fit_fun(ff, ud1, ud2);
			
			if (m2d(X_curr.y) <= m2d(X_next.y))
			{
				if (d > 0)
				{
					p[0] = m2d(X_prev.x);
					p[1] = m2d(X_next.x);
				}
				else
				{
					p[0] = m2d(X_next.x);
					p[1] = m2d(X_prev.x);
				}
				break;
			}
			X_prev = X_curr;
			X_curr = X_next;
			i++;
		}
		return p;
	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		int k = 1;
		double sqrt5 = sqrt(5.0);
		
		while (true)
		{
			double Fk = (pow((1 + sqrt5) / 2, k) - pow((1 - sqrt5) / 2, k)) / sqrt5;
			if (Fk > (b - a) / epsilon)
				break;
			k++;
		}
		
		int n = k;
		double a_i = a, b_i = b;
		
		double* F = new double[n + 1];
		F[0] = 1; F[1] = 1;
		for (int i = 2; i <= n; ++i)
			F[i] = F[i - 1] + F[i - 2];
			
		double c = b_i - (F[n - 1] / F[n]) * (b_i - a_i);
		double d = a_i + b_i - c;
		
		solution Xc, Xd;
		Xc.x = c; Xc.fit_fun(ff, ud1, ud2);
		Xd.x = d; Xd.fit_fun(ff, ud1, ud2);
		
		for (int i = 0; i < n - 3; ++i)
		{
			if (m2d(Xc.y) < m2d(Xd.y))
			{
				b_i = d;
				d = c;
				Xd = Xc;
				c = b_i - (F[n - i - 2] / F[n - i - 1]) * (b_i - a_i);
				Xc.x = c; Xc.fit_fun(ff, ud1, ud2);
			}
			else
			{
				a_i = c;
				c = d;
				Xc = Xd;
				d = a_i + b_i - c;
				Xd.x = d; Xd.fit_fun(ff, ud1, ud2);
			}
		}
		delete[] F;
		
		Xopt = Xc;
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}

}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		double c = (a + b) / 2.0;
		Xopt.x = a;
		Xopt.fit_fun(ff, ud1, ud2);
		double fa = m2d(Xopt.y);
		
		Xopt.x = b;
		Xopt.fit_fun(ff, ud1, ud2);
		double fb = m2d(Xopt.y);
		
		Xopt.x = c;
		Xopt.fit_fun(ff, ud1, ud2);
		double fc = m2d(Xopt.y);
		
		double d = c, d_prev = c, fd = fc;
		int i = 0;
		
		Xopt.ud = b - a;
		
		while (true) 
		{
			double l = fa * (b * b - c * c) + fb * (c * c - a * a) + fc * (a * a - b * b);
			double m = fa * (b - c) + fb * (c - a) + fc * (a - b);
			
			if (m <= 0.0) 
				throw string("m <= 0");
				
			d_prev = d;
			d = 0.5 * l / m;
			
			if (d <= a || d >= b || d == c) 
				throw string("d poza przedzialem");
				
			Xopt.x = d;
			Xopt.fit_fun(ff, ud1, ud2);
			fd = m2d(Xopt.y);
			
			if (d < c) 
			{
				if (fd < fc) 
				{
					b = c;
					fb = fc;
					c = d;
					fc = fd;
				}
				else 
				{
					a = d;
					fa = fd;
				}
			}
			else 
			{
				if (fd < fc) 
				{
					a = c;
					fa = fc;
					c = d;
					fc = fd;
				}
				else 
				{
					b = d;
					fb = fd;
				}
			}
			
			i++;
			Xopt.ud.add_row(b - a);
			
			if (solution::f_calls > Nmax) 
				throw string("Przekroczono Nmax");
				
			if ((b - a) < epsilon || (i > 1 && fabs(d - d_prev) < gamma)) 
				break;
		}
		
		Xopt.x = d;
		Xopt.y = fd;

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;                    // Obiekt przechowujący wynik optymalizacji

        Xopt.ud = trans(x0);             // Zapisujemy punkt początkowy do historii punktów

        solution XB(x0), XB_old, X;
        // XB    - aktualny punkt bazowy (najlepszy znaleziony punkt)
        // XB_old- poprzedni punkt bazowy
        // X     - punkt próbny, który będziemy sprawdzać

        XB.fit_fun(ff, ud1, ud2);         // Obliczamy wartość funkcji celu w punkcie początkowym

        while (true)                     // Główna pętla algorytmu
        {
            // ETAP 1: POSZUKIWANIE PRÓBNE
			// Sprawdzamy, czy można poprawić wynik, przesuwaąc się
    		// o krok s wzdłuż poszczególnych osi układu współrzędnych.
    		// Robi to funkcja HJ_trial().

    		X = HJ_trial(ff, XB, s, ud1, ud2);

    		if (X.y < XB.y)              // Jeżeli znaleziono punkt lepszy od bazowego
    		{
				// ETAP 2: RUCH ROBOCZY
        	// Skoro znaleźliśmy poprawę, próbujemy kontynuować
        	// poszukiwanie w kierunku, w którym funkcja malała.

        		while (true)                    // Zapamiętujemy poprzedni punkt bazowy
        		{                               // Nowy, lepszy punkt staje się punktem bazowym
            		XB_old = XB;
            		XB = X;

            		// Zapisujemy współrzędne nowego punktu bazowego do historii
            		// trans() zamienia wektor kolumnowy na wiersz,
            		// dzięki czemu każdy punkt zajmuje osobny wiersz macierzy.
            		Xopt.ud.add_row(trans(XB.x));

            		// Wyznaczamy nowy punkt w kierunku poprzedniej poprawy:
            		//
            		// X = XB + (XB - XB_old)
            		// X = 2*XB - XB_old
            		//
            		// Przykład:
            		// XB_old = (1, 2)
            		// XB     = (2, 3)
            		// X      = (3, 4)
            		//
					// Przesuwamy się dalej w tym kierunku
       				X.x = 2.0 * XB.x - XB_old.x;

        			X.fit_fun(ff, ud1, ud2); // Obliczamy funkcję celu w nowym punkcie

        			// Wykonujemy kolejne poszukiwanie próbne wokół nowego punktu.
        			// Może ono jeszcze bardziej poprawić znalezione rozwiązanie.
        			X = HJ_trial(ff, X, s, ud1, ud2);

        			// Jeżeli nowy punkt nie jest lepszy od aktualnego bazowego,
        			// kończymy ruch roboczy i wracamy do głównej pętli.
        			if (X.y >= XB.y)
        	    		break;

        			// Sprawdzamy, czy przekroczono dopuszczalną
        			// liczbę wywołań funkcji celu.
        			if (solution::f_calls > Nmax)
					{
        			Xopt = XB;                     // Zwracamy najlepszy znaleziony punkt
        			Xopt.flag = 0;                 // Flaga 0 oznacza przekroczenie limitu
        			return Xopt;                   // Kończymy optymalizację
    				}
				}
    		
			}
			else
    			// Jeżeli poszukiwanie próbne nie przyniosło poprawy,
    			// zmniejszamy długość kroku.
    			//
    			// Przykład:
				// s = 0.5, alpha = 0.5
        		// nowy krok: s = 0.25
        		// 
        		// Dzięki temu możemy dokładniej przeszukiwać
        		// otoczenie aktualnego punktu bazowego.
        		s *= alpha;

        	// PIERWSZE KRYTERIUM STOPU
        	// Jeżeli krok jest mniejszy od zadanej dokładności,
        	// uznajemy, że znaleźliśmy rozwiązanie.
        	if (s < epsilon)
        	{
    	    	Xopt = XB;                  // Zapisujemy najlepszy znaleziony punkt
    	    	Xopt.flag = 1;              // Flaga 1 oznacza prawidłowe zakończenie
            	break;
    		}

			// DRUGIE KRYTERIUM STOPU
    		// Kończymy, jeżeli wykorzystaliśmy limit wywołań funkcji celu
    		if (solution::f_calls > Nmax)
        	{
        		Xopt = XB;
            	Xopt.flag = 0;              // Nie osiągnięto kryterium dokładności
    	    	break;
    		}
		}

    	return Xopt;                        // Zwracamy wynik optymalizacji

	}
	catch (string ex_info)
	{
		// Przechwytujemy błąd i dodajemy informację,
        // że wystąpił on podczas działania funkcji HJ().
		throw ("solution HJ(...):\n" + ex_info);
	}
}

// Funkcja pomocnicza metody Hooke'a-Jeevesa
// HJ_trial() realizuje POSZUKIWANIE PRÓBNE.
//
// XB - punkt bazowy, od którego rozpoczynamy poszukiwanie
// s - długość kroku
//
// Dla każdej zmiennej sprawdzamy najpierw kierunek dodatni,
// a jeżeli nie ma poprawy, kierunek ujemny.
//
// Dla funkcji dwóch zmiennych sprawdzamy kolejno:
// (x1+s, x2), (x1-s, x2), (x1, x2+s), (x1, x2-s).
//
// Nie zawsze sprawdzane są wszystkie cztery punkty,
// ponieważ po znalezieniu poprawy nie sprawdzamy
// przeciwnego kierunku dla tej samej współrzędnej.
solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		int n = get_dim(XB);                    // Liczba zmiennych funkcji celu

        matrix D = ident_mat(n);                // Macierz jednostkowa n x n
        // Jej kolumny określają podstawowe kierunki poszukiwania.
		// 
        // Dla dwóch zmiennych:
        // 
        //       [1 0]
        // D = [ 0 1 ]
        // 
        // D[0] = (1, 0) - kierunek pierwszej współrzędnej
        // D[1] = (0, 1) - kierunek drugiej współrzędnej

        solution X;                    // Punkt próbny

        // Sprawdzamy kolejno każdy kierunek poszukiwania.
        for (int i = 0; i < n; ++i)
        {
            // KROK W KIERUNKU DODATNIM
            //
			// 
        	// Dla i = 0:
        	// X = XB + s*(1,0)
        	// 
        	// Dla i = 1:
        	// X = XB + s*(0,1)

        	X.x = XB.x + s * D[i];

        	X.fit_fun(ff, ud1, ud2); // Obliczamy wartość funkcji w punkcie próbnym

        	// Jeżeli wartość funkcji jest mniejsza,
        	// przyjmujemy punkt próbny jako nowy punkt bazowy.
        	if (X.y < XB.y)
        	    XB = X;
        	else
        	{
        	    // KROK W KIERUNKU UJEMNYM
        	    //
				// Jeżeli ruch dodatni nie poprawił wyniku,
            	// sprawdzamy przeciwny kierunek.

        	    X.x = XB.x - s * D[i];

        	    X.fit_fun(ff, ud1, ud2);

        	    // Jeżeli ruch ujemny poprawił wynik,
        	    // aktualizujemy punkt bazowy.
        	    if (X.y < XB.y)
        	        XB = X;
        	}

        	// Po sprawdzeniu jednego kierunku przechodzimy
        	// do kolejnej zmiennej.
        	// 
        	// Ważne: jeżeli znaleźliśmy poprawę,
        	// kolejne kroki wykonujemy już od nowego XB.
		}
		// Zwracamy najlepszy punkt znaleziony
        // podczas poszukiwania próbnego.
		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

// Metoda Rosenbrocka (Rosen)
// Jest to metoda bezgradientowa, która oprócz zmiany
// długości kroków potrafi również zmieniać kierunki poszukiwania.
// 
// W metodzie Hooke'a-Jeevesa poszukujemy wzdłuż stałych osi.
// W metodzie Rosenbrocka kierunki mogą zostać obrócone,
// aby lepiej odpowiadały kształtowi optymalizowanej funkcji.
// 
// ff        - funkcja celu
// x0        - punkt początkowy
// s0        - wektor początkowych długości kroków
// alpha     - współczynnik zwiększania kroku po poprawie (np. 2)
// beta      - współczynnik zmniejszania kroku po niepowodzeniu (np. 0.5)
// epsilon   - zadana dokładność
// Nmax      - maksymalna liczba wywołań funkcji celu
// ud1, ud2  - dodatkowe dane funkcji celu
solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;                    // Obiekt przechowujący wynik

        Xopt.ud = trans(x0);             // Zapisujemy punkt początkowy do historii

		solution XB(x0), X;
        // XB - aktualny najlepszy punkt
        // X - punkt próbny

        int n = get_dim(XB);                    // Liczba zmiennych funkcji celu

        matrix l(n, 1), p(n, 1), s(s0), D = ident_mat(n);
        // l - suma udanych przesunięć w poszczególnych kierunkach
        // p - licznik nieudanych prób w poszczególnych kierunkach
        // s - aktualne długości kroków (osobny krok dla każdej zmiennej)
        // D - macierz kierunków poszukiwania
        // 
        // Początkowo D jest macierzą jednostkową,
        // czyli poszukujemy wzdłuż osi układu współrzędnych.
        // 
        // W odróżnieniu od HJ, gdzie mamy jeden krok s,
        // tutaj każdy kierunek ma własny krok s(i).

        XB.fit_fun(ff, ud1, ud2);               // Obliczamy funkcję celu w punkcie początkowym

		// Główna pętla metody Rosenbrocka
        while (true)
        {
            // ETAP 1: POSZUKIWANIE W AKTUALNYCH KIERUNKACH
            // 
            // Sprawdzamy kolejno każdy z n kierunków.

            for (int i = 0; i < n; ++i)
            {
                // Wyznaczamy nowy punkt:
                // 
                // X = XB + s(i)*D[i]
                // 
                // D[i] określa kierunek,
                // a s(i) określa długość i zwrot kroku.
                // 
                // Ujemne s(i) oznacza ruch w przeciwną stronę.

                X.x = XB.x + s(i) * D[i];

                X.fit_fun(ff, ud1, ud2); // Obliczamy funkcję celu

				if (X.y < XB.y)
                {
                    // SUKCES - znaleziono lepszy punkt
                    // Aktualizujemy najlepszy punkt

                    XB = X;

                    // Sumujemy udane przesunięcia
                    l(i) += s(i);
                    // w kierunku i

                    // Zwiększamy długość kroku
                    s(i) *= alpha;
					// w tym kierunku
// 
// Przykład:
// s(i) = 0.5, alpha = 2
// nowy krok: s(i) = 1.0
//  
// Skoro ruch poprawił wynik,
// próbujemy wykonywać większe kroki.              
				}
                else
                {
                    // NIEPOWODZENIE - nie znaleziono lepszego punktu

                    ++p(i);                  // Zwiększamy licznik niepowodzeń

                    s(i) *= -beta;           // Zmniejszamy długość kroku
					// i zmieniamy jego znak
// 
// Przykład:
// s(i) = 1.0, beta = 0.5
// nowy krok: s(i) = -0.5
// 
// W kolejnej próbie sprawdzimy przeciwny kierunek,
// wykonując jednocześnie krótszy krok.
                }
            }	
			
			// Zapisujemy aktualny najlepszy punkt do historii
        	// Dla funkcji dwóch zmiennych każdy wiersz
        	// zawiera współrzędne (x1, x2).
        	Xopt.ud.add_row(trans(XB.x));

        	// ETAP 2: SPRAWDZENIE MOŻLIWOŚCI ZMIANY KIERUNKÓW
        	// 
        	// Nowe kierunki wyznaczamy dopiero wtedy,
        	// gdy w KAŻDYM kierunku wystąpiły:
        	// 
        	// - udane przesunięcia (l(i) != 0),
        	// - nieudane próby (p(i) != 0).

        	bool change = true;

        	for (int i = 0; i < n; ++i)
        	    if (l(i) == 0 || p(i) == 0)
				{
					// Jeżeli w którymś kierunku nie spełniono
            		// obu warunków, nie zmieniamy jeszcze kierunków.
            		change = false;
            		break;
				}
			
			if (change)
        	{
        	    // ETAP 3: WYZNACZANIE NOWYCH KIERUNKÓW
            	// 
            	// Budujemy nowy układ kierunków poszukiwania.
            	// Wykorzystujemy informacje o udanych przesunięciach,
            	// zapisane w wektorze l.

            	matrix Q(n, n), v(n, 1);

				// Tworzymy pomocniczą macierz Q.
				// Wypełniamy jej dolną część wartościami l(i).
				for (int i = 0; i < n; ++i)
					for (int j = 0; j <= i; ++j)
						Q(i, j) = l(i);

				// Przekształcamy kierunki za pomocą aktualnej macierzy D.
				// Dzięki temu nowe kierunki uwzględniają dotychczasowy
				// układ współrzędnych poszukiwania.

				Q = D * Q;

				// Pierwszy nowy kierunek otrzymujemy przez
				// normalizację pierwszego wektora Q.
				// 
				// norm() oblicza długość wektora.
				// Dzielenie przez długość daje wektor jednostkowy.

				v = Q[0] / norm(Q[0]);

				// Zapisujemy pierwszy nowy kierunek w macierzy D.
				D.set_col(v, 0);

				// Kolejne kierunki wyznaczamy metodą
				// ortogonalizacji Grama-Schmidta.
				// 
				// Jej zadaniem jest utworzenie wektorów:
				// - prostopadłych do siebie,
				// - o długości równej 1.
				// 
				// Dzięki temu otrzymujemy nowy układ kierunków
				// dostosowany do przebiegu dotychczasowych poszukiwań.

				for (int i = 1; i < n; ++i)
				{
					matrix temp(n, 1);

					// Obliczamy sumę rzutów aktualnego wektora Q[i]
					// na wcześniej wyznaczone kierunki D[j].
					// 
					// Odejmując tę sumę, usuwamy składowe
					// równoległe do poprzednich kierunków.

					for (int j = 0; j < i; ++j)
						temp = temp + trans(Q[i]) * D[j] * D[j];

					// Wyznaczamy składową prostopadłą do poprzednich kierunków
					// i normalizujemy ją do długości 1.

					v = (Q[i] - temp) / norm(Q[i] - temp);

					// Zapisujemy nowy kierunek w odpowiedniej kolumnie D.
					D.set_col(v, i);
				}

				// Po zmianie kierunków rozpoczynamy kolejny etap poszukiwań.
				// 
				// Przywracamy początkowe długości kroków.
				s = s0;

				// Zerujemy informacje o udanych przesunięciach.
				l = matrix(n, 1);

				// Zerujemy liczniki niepowodzeń.
				p = matrix(n, 1);
			}

			// ETAP 4: SPRAWDZENIE KRYTERIÓW STOPU

			// Szukamy największej bezwzględnej długości kroku.
			// 
			// Używamy abs(), ponieważ krok może być ujemny.
			// Interesuje nas jego długość, a nie zwrot.

			double max_s = abs(s(0));

			for (int i = 1; i < n; ++i)
				if (max_s < abs(s(i)))
					max_s = abs(s(i));

			// PIERWSZE KRYTERIUM STOPU
			// Jeżeli wszystkie kroki są mniejsze od epsilon,
			// kończymy optymalizację.
			// 
			// Wystarczy sprawdzić największy krok.
			// Jeżeli największy jest mniejszy od epsilon,
			// to pozostałe również.

			if (max_s < epsilon)
			{
				Xopt = XB;                  // Zwracamy najlepszy znaleziony punkt
				Xopt.flag = 1;              // Flaga 1 - prawidłowe zakończenie
				break;
			}

			// DRUGIE KRYTERIUM STOPU
			// Kończymy, gdy przekroczymy limit
			// wywołań funkcji celu.

			if (solution::f_calls > Nmax)
			{
				Xopt = XB;
				Xopt.flag = 0;              // Flaga 0 - przekroczony limit
				break;
			}
		}

		return Xopt;						// Zwracamy wynik optymalizacji
	}
	catch (string ex_info)
	{
		// Obsługa błędów występujących podczas działania Rosen().
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
