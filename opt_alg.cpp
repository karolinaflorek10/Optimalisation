
#include "opt_alg.h"
#include <vector>
#include <cmath>
#include <algorithm>

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
    // Zmienne wejsciowe:
    // ff - wskaznik do funkcji celu
    // N - liczba zmiennych funkcji celu
    // lb, ub - dolne i gorne ograniczenie
    // epsilon - zakladana dokladnosc rozwiazania
    // Nmax - maksymalna liczba wywolan funkcji celu
    // ud1, ud2 - user data
    try
    {
       solution Xopt;
       Xopt = rand_mat(N);
       for (int i = 0; i < N; ++i)
          Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
       Xopt.fit_fun(ff, ud1, ud2);
       Xopt.ud = (Xopt.y);
       if (Xopt.y < epsilon)
       {
          Xopt.flag = 1;
          return Xopt;
       }
       solution X;
       while (true)
       {
          X = rand_mat(N);
          for (int i = 0; i < N; ++i)
             X.x(i) = (ub(i) - lb(i)) * X.x(i) + lb(i);
          X.fit_fun(ff, ud1, ud2);
          if (X.y < Xopt.y)
          {
             Xopt = X;
             if (Xopt.y < epsilon)
             {
                Xopt.flag = 1;
                break;
             }
          }
          if (solution::f_calls > Nmax)
          {
             Xopt.flag = 0;
             break;
          }
          Xopt.ud.add_row(Xopt.y);
       }
       Xopt.ud.add_row(Xopt.y);
       return Xopt;
    }
    catch (string ex_info)
    {
       throw ("solution MC(...):\n" + ex_info);
    }
}


// ################## LAB 1 #################
// Metoda ekspansji
double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
    try
    {
       double* p = new double[2] { 0, 0 };
       solution x_prev(x0), x_curr(x0 + d);

       x_prev.fit_fun(ff, ud1, ud2);
       x_curr.fit_fun(ff, ud1, ud2);

       if (x_curr.y == x_prev.y)
       {
          p[0] = x_prev.x(0);
          p[1] = x_curr.x(0);
          if (p[0] > p[1]) swap(p[0], p[1]);
          return p;
       }

       if (x_curr.y > x_prev.y)
       {
          d = -d;
          x_curr.x = x0 + d;
          x_curr.fit_fun(ff, ud1, ud2);

          if (x_curr.y >= x_prev.y)
          {
             p[0] = x_curr.x(0);
             p[1] = x0 - d;
             if (p[0] > p[1]) swap(p[0], p[1]);
             return p;
          }
       }

       solution x_next;
       int i = 0;
       while (true)
       {
          if (solution::f_calls > Nmax)
             throw string("expansion: przekroczono maksymalna liczbe wywolan Nmax");

          i++;
          x_next.x = x0 + pow(alpha, i) * d;
          x_next.fit_fun(ff, ud1, ud2);

          if (x_curr.y <= x_next.y)
             break;

          x_prev = x_curr;
          x_curr = x_next;
       }

       if (d > 0)
       {
          p[0] = x_prev.x(0);
          p[1] = x_next.x(0);
       }
       else
       {
          p[0] = x_next.x(0);
          p[1] = x_prev.x(0);
       }

       if (p[0] > p[1]) swap(p[0], p[1]);
       return p;
    }
    catch (string ex_info)
    {
       throw ("double* expansion(...):\n" + ex_info);
    }
}


// Metoda Fibonacciego
solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
    try
    {
       solution Xopt;
       bool pokaz_zbieznosc = (a == -100.0 && b == 100.0);

       if (pokaz_zbieznosc)
       {
          cout << "\nFIBONACCI ZBIEZNOSC\n";
          cout << "0 " << b - a << endl;
       }

       vector<double> F = { 1.0, 1.0 };

       while (F.back() <= (b - a) / epsilon)
          F.push_back(F.back() + F[F.size() - 2]);

       int k = static_cast<int>(F.size()) - 1;

       solution c, d;
       c.x = b - (F[k - 1] / F[k]) * (b - a);
       d.x = a + b - c.x(0);

       c.fit_fun(ff, ud1, ud2);
       d.fit_fun(ff, ud1, ud2);

       for (int i = 0; i <= k - 3; ++i)
       {
          if (c.y < d.y)
          {
             b = d.x(0);
             d = c;
             c.x = b - (F[k - i - 2] / F[k - i - 1]) * (b - a);
             c.fit_fun(ff, ud1, ud2);
          }
          else
          {
             a = c.x(0);
             c = d;
             d.x = a + b - c.x(0);
             d.fit_fun(ff, ud1, ud2);
          }
          if (pokaz_zbieznosc)
              cout << i + 1 << " " << b - a << endl;
       }

       Xopt = (c.y < d.y) ? c : d;
       Xopt.flag = 1;
       return Xopt;
    }
    catch (string ex_info)
    {
       throw ("solution fib(...):\n" + ex_info);
    }
}


// Metoda Lagrange'a
solution lag(matrix(*ff)(matrix, matrix, matrix),
    double a, double b, double epsilon, double gamma,
    int Nmax, matrix ud1, matrix ud2)
{
    try
    {
        if (a >= b || epsilon <= 0 || gamma <= 0 || Nmax <= 0)
            throw string("lag: niepoprawne parametry");

        solution A(a), B(b), C((a + b) / 2.0), D;

        A.fit_fun(ff, ud1, ud2);
        B.fit_fun(ff, ud1, ud2);
        C.fit_fun(ff, ud1, ud2);

        double poprzednie_d = NAN;
        bool pokaz_zbieznosc = (a == -100.0 && b == 100.0);
        int iteracja = 0;

        if (pokaz_zbieznosc)
        {
            cout << "\nLAGRANGE ZBIEZNOSC\n";
            cout << "0 " << b - a << endl;
        }

        while (true)
        {
            if (solution::f_calls >= Nmax)
                throw string("lag: przekroczono Nmax");

            double av = A.x(0);
            double bv = B.x(0);
            double cv = C.x(0);

            // Wspolczynniki interpolacji
            double l =
                A.y(0) * (bv * bv - cv * cv) +
                B.y(0) * (cv * cv - av * av) +
                C.y(0) * (av * av - bv * bv);

            double m =
                A.y(0) * (bv - cv) +
                B.y(0) * (cv - av) +
                C.y(0) * (av - bv);

            if (m <= 0 || !std::isfinite(m))
                throw string("lag: niepoprawna interpolacja");

            double nowe_d = 0.5 * l / m;

            if (!std::isfinite(nowe_d) ||
                nowe_d <= av || nowe_d >= bv)
                throw string("lag: punkt poza przedzialem");

            D.x = nowe_d;
            D.fit_fun(ff, ud1, ud2);

            if (nowe_d < cv)
            {
                if (D.y < C.y)
                {
                    B = C;
                    C = D;
                }
                else
                {
                    A = D;
                }
            }
            else if (nowe_d > cv)
            {
                if (D.y < C.y)
                {
                    A = C;
                    C = D;
                }
                else
                {
                    B = D;
                }
            }
            else
            {
                C = D;
                ++iteracja;
                if (pokaz_zbieznosc)
                    cout << iteracja << " " << B.x(0) - A.x(0) << endl;
                break;
            }

            ++iteracja;

            if (pokaz_zbieznosc)
                cout << iteracja << " " << B.x(0) - A.x(0) << endl;

            if (B.x(0) - A.x(0) < epsilon)
                break;

            if (std::isfinite(poprzednie_d) &&
                std::abs(nowe_d - poprzednie_d) < gamma)
                break;

            poprzednie_d = nowe_d;
        }

        solution Xopt = (D.y < C.y) ? D : C;
        Xopt.flag = 1;

        return Xopt;
    }
    catch (string ex_info)
    {
        throw ("solution lag(...):\n" + ex_info);
    }
}


// ################## LAB 2 #################
// Metoda Hooke'a-Jeevesa
solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0,
    double s, double alpha, double epsilon, int Nmax,
    matrix ud1, matrix ud2)
{
    try
    {
        solution Xopt;
        Xopt.ud = trans(x0);

        solution XB(x0), XB_old, X;
        XB.fit_fun(ff, ud1, ud2);

        do
        {
            X = HJ_trial(ff, XB, s, ud1, ud2);

            if (X.y < XB.y)
            {
                do
                {
                    XB_old = XB;
                    XB = X;
                    Xopt.ud.add_row(trans(XB.x));

                    X.x = 2.0 * XB.x - XB_old.x;
                    X.fit_fun(ff, ud1, ud2);
                    X = HJ_trial(ff, X, s, ud1, ud2);

                    if (solution::f_calls > Nmax)
                    {
                        matrix historia = Xopt.ud;
                        Xopt = XB;
                        Xopt.ud = historia;
                        Xopt.flag = 0;
                        return Xopt;
                    }

                } while (X.y < XB.y);
            }
            else
            {
                s = alpha * s;
            }

            if (solution::f_calls > Nmax)
            {
                matrix historia = Xopt.ud;
                Xopt = XB;
                Xopt.ud = historia;
                Xopt.flag = 0;
                return Xopt;
            }

        } while (s >= epsilon);

        matrix historia = Xopt.ud;
        Xopt = XB;
        Xopt.ud = historia;
        Xopt.flag = 1;

        return Xopt;
    }
    catch (string ex_info)
    {
        throw ("solution HJ(...):\n" + ex_info);
    }
}


// Kroki probne metody Hooke'a-Jeevesa
solution HJ_trial(matrix(*ff)(matrix, matrix, matrix),
    solution XB, double s, matrix ud1, matrix ud2)
{
    try
    {
        int n = get_dim(XB);
        matrix E = ident_mat(n);
        solution X;

        for (int j = 0; j < n; ++j)
        {
            X.x = XB.x + s * E[j];
            X.fit_fun(ff, ud1, ud2);

            if (X.y < XB.y)
            {
                XB = X;
            }
            else
            {
                X.x = XB.x - s * E[j];
                X.fit_fun(ff, ud1, ud2);

                if (X.y < XB.y)
                {
                    XB = X;
                }
            }
        }

        return XB;
    }
    catch (string ex_info)
    {
        throw ("solution HJ_trial(...):\n" + ex_info);
    }
}


// Metoda Rosenbrocka
solution Rosen(matrix(*ff)(matrix, matrix, matrix),
    matrix x0, matrix s0, double alpha, double beta,
    double epsilon, int Nmax, matrix ud1, matrix ud2)
{
    try
    {
        solution Xopt;
        Xopt.ud = trans(x0);

        solution XB(x0), X;
        int n = get_dim(XB);

        matrix D = ident_mat(n);
        matrix lambda(n, 1);
        matrix p(n, 1);
        matrix s(s0);

        XB.fit_fun(ff, ud1, ud2);

        while (true)
        {
            for (int j = 0; j < n; ++j)
            {
                X.x = XB.x + s(j) * D[j];
                X.fit_fun(ff, ud1, ud2);

                if (X.y < XB.y)
                {
                    XB = X;
                    lambda(j) += s(j);
                    s(j) *= alpha;
                }
                else
                {
                    ++p(j);
                    s(j) *= -beta;
                }
            }

            Xopt.ud.add_row(trans(XB.x));

            bool change_basis = true;

            for (int j = 0; j < n; ++j)
            {
                if (lambda(j) == 0 || p(j) == 0)
                {
                    change_basis = false;
                    break;
                }
            }

            if (change_basis)
            {
                matrix Q(n, n), v(n, 1);

                for (int j = 0; j < n; ++j)
                {
                    for (int k = 0; k <= j; ++k)
                    {
                        Q(j, k) = lambda(j);
                    }
                }

                Q = D * Q;

                if (norm(Q[0]) > 1e-12)
                {
                    v = Q[0] / norm(Q[0]);
                    D.set_col(v, 0);

                    for (int j = 1; j < n; ++j)
                    {
                        matrix temp(n, 1);

                        for (int k = 0; k < j; ++k)
                        {
                            temp = temp + trans(Q[j]) * D[k] * D[k];
                        }

                        matrix pozostaly = Q[j] - temp;

                        if (norm(pozostaly) > 1e-12)
                        {
                            v = pozostaly / norm(pozostaly);
                            D.set_col(v, j);
                        }
                    }
                }

                lambda = matrix(n, 1);
                p = matrix(n, 1);
                s = s0;
            }

            if (solution::f_calls > Nmax)
            {
                matrix historia = Xopt.ud;
                Xopt = XB;
                Xopt.ud = historia;
                Xopt.flag = 0;
                return Xopt;
            }

            double max_s = abs(s(0));

            for (int j = 1; j < n; ++j)
            {
                if (abs(s(j)) > max_s)
                {
                    max_s = abs(s(j));
                }
            }

            if (max_s < epsilon)
            {
                matrix historia = Xopt.ud;
                Xopt = XB;
                Xopt.ud = historia;
                Xopt.flag = 1;
                return Xopt;
            }
        }
    }
    catch (string ex_info)
    {
        throw ("solution Rosen(...):\n" + ex_info);
    }
}


// ################## LAB 3 #################

// Metoda funkcji kary.
// Kazde minimum czastkowe liczy metoda Neldera-Meada.
solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0,
             double c, double dc, double epsilon, int Nmax,
             matrix ud1, matrix ud2)
{
    if (c <= 0.0 || dc <= 0.0 || epsilon <= 0.0 || Nmax <= 0)
        throw string("pen: niepoprawne parametry");

    const double alpha = 1.0;
    const double beta = 0.5;
    const double gamma = 2.0;
    const double delta = 0.5;
    const double s = 0.5;

    solution previous(x0), current;

    while (true)
    {
        if (solution::f_calls >= Nmax)
        {
            previous.flag = 0;
            return previous;
        }

        current = sym_NM(ff, previous.x, s, alpha, beta, gamma,
                         delta, epsilon, Nmax, ud1, c);

        // Nie ignorujemy niepowodzenia wewnetrznej optymalizacji.
        if (current.flag != 1)
        {
            current.flag = 0;
            return current;
        }

        if (norm(current.x - previous.x) < epsilon)
        {
            current.flag = 1;
            return current;
        }

        previous = current;
        c *= dc;

        if (!std::isfinite(c) || c <= 0.0)
            throw string("pen: wspolczynnik kary poza zakresem");
    }
}


// Metoda sympleks Neldera-Meada
solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0,
                double s, double alpha, double beta, double gamma,
                double delta, double epsilon, int Nmax,
                matrix ud1, matrix ud2)
{
    if (s <= 0 || epsilon <= 0 || Nmax <= 0)
        throw string("sym_NM: niepoprawne parametry");

    const int n = get_len(x0);
    const int N = n + 1;

    matrix D = ident_mat(n);
    std::vector<solution> S(N);

    auto evaluate = [&](solution& point)
    {
        point.fit_fun(ff, ud1, ud2);

        if (!std::isfinite(m2d(point.y)))
            point.y = 1e15;
    };

    // Tworzenie sympleksu
    S[0].x = x0;
    evaluate(S[0]);

    for (int i = 1; i < N; ++i)
    {
        S[i].x = x0 + s * D[i - 1];
        evaluate(S[i]);
    }

    while (true)
    {
        // Wyznaczamy najlepszy i najgorszy wierzcholek.
        // Indeksy sa rozne nawet przy remisie.
        std::vector<int> order(N);

        for (int i = 0; i < N; ++i)
            order[i] = i;

        std::sort(order.begin(), order.end(),
            [&](int a, int b)
            {
                return m2d(S[a].y) < m2d(S[b].y);
            });

        const int best = order.front();
        const int worst = order.back();

        if (solution::f_calls >= Nmax)
        {
            solution result = S[best];
            result.flag = 0;
            return result;
        }

        // Kryterium stopu
        double span = 0.0;

        for (int i = 0; i < N; ++i)
            span = std::max(span, norm(S[i].x - S[best].x));

        if (span < epsilon)
        {
            solution result = S[best];
            result.flag = 1;
            return result;
        }

        // Srodek ciezkosci bez najgorszego wierzcholka
        matrix centroid(n, 1);

        for (int i = 0; i < N; ++i)
            if (i != worst)
                centroid = centroid + S[i].x;

        centroid = centroid / double(n);

        // Odbicie
        solution reflection;
        reflection.x = centroid + alpha * (centroid - S[worst].x);
        evaluate(reflection);

        if (reflection.y < S[best].y)
        {
            // Ekspansja
            solution expansion;
            expansion.x = centroid + gamma * (reflection.x - centroid);
            evaluate(expansion);

            S[worst] = (expansion.y < reflection.y)
                         ? expansion : reflection;
        }
        else if (reflection.y < S[worst].y)
        {
            S[worst] = reflection;
        }
        else
        {
            // Zawezenie
            solution contraction;
            contraction.x = centroid + beta * (S[worst].x - centroid);
            evaluate(contraction);

            if (contraction.y < S[worst].y)
            {
                S[worst] = contraction;
            }
            else
            {
                // Redukcja sympleksu
                for (int i = 0; i < N; ++i)
                {
                    if (i == best)
                        continue;

                    S[i].x = S[best].x +
                             delta * (S[i].x - S[best].x);

                    evaluate(S[i]);
                }
            }
        }

        // Po zmianie sympleksu indeksy zostana
        // wyznaczone ponownie w nastepnej iteracji.
    }
}






// ################## KOLEJNE LABORATORIA ##################

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


// ----------------- LAB 4 ------------------------
//--------------------------------------------------
//metoda najszybszego spadku
solution SD(matrix(*ff)(matrix, matrix, matrix),
	matrix(*gf)(matrix, matrix, matrix),
	matrix x0, double h0, double epsilon, int Nmax,
	matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;

		// Zapisujemy punkt poczatkowy do historii.
		// Dane z pola ud beda wykorzystywane do wykresow.
		Xopt.ud = trans(x0);

		// X0 - aktualny punkt
		// X1 - nowy punkt
		// h  - optymalna dlugosc kroku
		solution X0(x0), X1, h;

		// Liczba zmiennych funkcji celu.
		int n = get_len(x0);

		// d - kierunek poszukiwania
		// P - macierz przechowujaca punkt i kierunek
		matrix d(n, 1), P(n, 2);

		double* ab = nullptr;

		while (true)
		{
			// KROK 1:
			// Wyznaczamy kierunek przeciwny do gradientu.
			d = -X0.grad(gf, ud1, ud2);

			// KROK 2:
			// Wyznaczamy dlugosc kroku.
			if (h0 < 0)
			{
				// WERSJA ZMIENNOKROKOWA

				// Zapisujemy punkt poczatkowy
				// i kierunek poszukiwania.
				P.set_col(X0.x, 0);
				P.set_col(d, 1);

				// Najpierw znajdujemy przedzial,
				// w ktorym znajduje sie minimum.
				ab = expansion(
					ff, 0, 1, 1.2,
					Nmax, ud1, P
				);

				// Nastepnie metoda zlotego podzialu
				// znajduje najlepsza dlugosc kroku.
				h = golden(
					ff, ab[0], ab[1],
					epsilon, Nmax, ud1, P
				);

				// Zwalniamy pamiec po expansion().
				delete[] ab;
				ab = nullptr;

				// Aktualizujemy punkt.
				X1.x = X0.x + h.x * d;
			}
			else
			{
				// WERSJA STALOKROKOWA

				// Korzystamy z kroku podanego
				// przez uzytkownika, np. 0.05 lub 0.25.
				X1.x = X0.x + h0 * d;
			}

			// Zapisujemy nowy punkt do historii.
			Xopt.ud.add_row(trans(X1.x));

			// KROK 3:
			// Sprawdzamy kryterium stopu.
			if (norm(X1.x - X0.x) < epsilon)
			{
				Xopt = X1;

				// Obliczamy wartosc funkcji celu
				// w otrzymanym punkcie.
				Xopt.fit_fun(ff, ud1, ud2);

				// Poprawne zakonczenie optymalizacji.
				Xopt.flag = 1;
				break;
			}

			// KROK 4:
			// Sprawdzamy limit wywolan funkcji celu
			// i gradientu.
			if (solution::f_calls > Nmax ||
				solution::g_calls > Nmax)
			{
				Xopt = X1;
				Xopt.fit_fun(ff, ud1, ud2);

				// Przekroczono limit.
				Xopt.flag = 0;
				break;
			}

			// Przygotowujemy kolejna iteracje.
			X0 = X1;
		}

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

// metoda gradientów spężonych
solution CG(matrix(*ff)(matrix, matrix, matrix),
	matrix(*gf)(matrix, matrix, matrix),
	matrix x0, double h0, double epsilon, int Nmax,
	matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;

		// Historia kolejnych punktow.
		Xopt.ud = trans(x0);

		solution X0(x0), X1, h;
		int n = get_len(x0);
		matrix d(n, 1), P(n, 2);
		double* ab = nullptr;

		// Pierwszy kierunek jest taki sam
		// jak w metodzie najszybszego spadku.
		d = -X0.grad(gf, ud1, ud2);

		while (true)
		{
			// Wyznaczamy dlugosc kroku.
			if (h0 < 0)
			{
				// WERSJA ZMIENNOKROKOWA

				P.set_col(X0.x, 0);
				P.set_col(d, 1);

				// Znajdujemy przedzial poszukiwan.
				ab = expansion(
					ff, 0, 1, 1.2,
					Nmax, ud1, P
				);

				// Znajdujemy optymalna dlugosc kroku.
				h = golden(
					ff, ab[0], ab[1],
					epsilon, Nmax, ud1, P
				);

				delete[] ab;
				ab = nullptr;

				X1.x = X0.x + h.x * d;
			}
			else
			{
				// WERSJA STALOKROKOWA
				X1.x = X0.x + h0 * d;
			}

			// Zapisujemy nowy punkt.
			Xopt.ud.add_row(trans(X1.x));

			// KRYTERIUM STOPU
			if (norm(X1.x - X0.x) < epsilon)
			{
				Xopt = X1;
				Xopt.fit_fun(ff, ud1, ud2);

				Xopt.flag = 1;
				break;
			}

			// Sprawdzamy limit wywolan.
			if (solution::f_calls > Nmax ||
				solution::g_calls > Nmax)
			{
				Xopt = X1;
				Xopt.fit_fun(ff, ud1, ud2);

				Xopt.flag = 0;
				break;
			}

			// Wyznaczamy gradient w nowym punkcie.
			matrix g1 = X1.grad(gf, ud1, ud2);

			// Norma poprzedniego gradientu.
			double old_norm = norm(X0.g);

			// Zabezpieczenie przed dzieleniem przez zero.
			if (old_norm == 0.0)
			{
				Xopt = X1;
				Xopt.fit_fun(ff, ud1, ud2);

				Xopt.flag = 1;
				break;
			}

			// WSPOLCZYNNIK FLETCHERA-REEVESA
			double beta =
				pow(norm(g1), 2) /
				pow(old_norm, 2);

			// NOWY KIERUNEK SPRZEZONY
			d = -g1 + beta * d;

			// Aktualizujemy punkt.
			X0 = X1;
		}

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

//metoda Newtona

solution Newton(matrix(*ff)(matrix, matrix, matrix),
	matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix),
	matrix x0, double h0,
	double epsilon, int Nmax,
	matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;

		// Zapisujemy punkt poczatkowy.
		Xopt.ud = trans(x0);

		solution X0(x0), X1, h;
		int n = get_len(x0);
		matrix d(n, 1), P(n, 2);
		double* ab = nullptr;

		while (true)
		{
			// KROK 1:
			// Obliczamy hesjan i gradient.

			// Kierunek Newtona:
			// d = -H^(-1) * gradient
			d = -inv(X0.hess(Hf, ud1, ud2)) *
				X0.grad(gf, ud1, ud2);

			// KROK 2:
			// Wyznaczamy dlugosc kroku.
			if (h0 < 0)
			{
				// WERSJA ZMIENNOKROKOWA

				P.set_col(X0.x, 0);
				P.set_col(d, 1);

				ab = expansion(
					ff, 0, 1, 1.2,
					Nmax, ud1, P
				);

				h = golden(
					ff, ab[0], ab[1],
					epsilon, Nmax, ud1, P
				);

				delete[] ab;
				ab = nullptr;

				X1.x = X0.x + h.x * d;
			}
			else
			{
				// WERSJA STALOKROKOWA
				X1.x = X0.x + h0 * d;
			}

			// Zapisujemy historie punktow.
			Xopt.ud.add_row(trans(X1.x));

			// KRYTERIUM STOPU
			if (norm(X1.x - X0.x) < epsilon)
			{
				Xopt = X1;
				Xopt.fit_fun(ff, ud1, ud2);

				Xopt.flag = 1;
				break;
			}

			// LIMIT WYWOLAN
			if (solution::f_calls > Nmax ||
				solution::g_calls > Nmax ||
				solution::H_calls > Nmax)
			{
				Xopt = X1;
				Xopt.fit_fun(ff, ud1, ud2);

				Xopt.flag = 0;
				break;
			}

			// Kolejna iteracja.
			X0 = X1;
		}

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

//metoda złotego podziału

solution golden(matrix(*ff)(matrix, matrix, matrix),
	double a, double b,
	double epsilon, int Nmax,
	matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;

		// Wspolczynnik zlotego podzialu.
		double alfa = (sqrt(5.0) - 1.0) / 2.0;

		// A, B - granice przedzialu
		// C, D - punkty wewnetrzne
		solution A(a), B(b), C, D;

		// Wyznaczamy punkty wewnetrzne.
		C.x = B.x - alfa * (B.x - A.x);
		D.x = A.x + alfa * (B.x - A.x);

		// Obliczamy wartosci funkcji.
		C.fit_fun(ff, ud1, ud2);
		D.fit_fun(ff, ud1, ud2);

		while (true)
		{
			if (C.y < D.y)
			{
				// Minimum znajduje sie
				// w lewej czesci przedzialu.

				B = D;
				D = C;

				C.x = B.x - alfa * (B.x - A.x);

				C.fit_fun(ff, ud1, ud2);
			}
			else
			{
				// Minimum znajduje sie
				// w prawej czesci przedzialu.

				A = C;
				C = D;

				D.x = A.x + alfa * (B.x - A.x);

				D.fit_fun(ff, ud1, ud2);
			}

			// Sprawdzamy szerokosc przedzialu.
			if (B.x - A.x < epsilon)
			{
				// Zwracamy srodek koncowego przedzialu.
				Xopt.x = (A.x + B.x) / 2.0;

				Xopt.fit_fun(ff, ud1, ud2);

				Xopt.flag = 1;
				break;
			}

			// Sprawdzamy limit wywolan.
			if (solution::f_calls > Nmax)
			{
				Xopt.x = (A.x + B.x) / 2.0;

				Xopt.fit_fun(ff, ud1, ud2);

				Xopt.flag = 0;
				break;
			}
		}

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}
