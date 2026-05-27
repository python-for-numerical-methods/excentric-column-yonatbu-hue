def find_critical_load(L, E, A, r, c, e, sigma_allow):
    """
    L: אורך במ"מ
    E: מודול אלסטיות ב-MPa
    A: שטח חתך בממ"ר
    r: רדיוס אינרציה במ"מ
    c: מרחק לסיב קיצוני במ"מ
    e: אקסצנטריות במ"מ
    sigma_allow: מאמץ מותר ב-MPa

    Return: העומס P בניוטון (float)
    """
 def sigma_max(P):
        """
        Maximum stress according to the secant formula.
        """

        # Argument of secant function
        theta = (L / (2 * r)) * np.sqrt(P / (E * A))

        # sec(x) = 1 / cos(x)
        sec_theta = 1 / np.cos(theta)

        return (P / A) * (1 + (e * c / r**2) * sec_theta)

    def f(P):
        """
        Root function:
        sigma_max(P) - allowable stress
        """
        return sigma_max(P) - sigma_allow

    # Euler buckling load (used as upper bound)
    I = A * r**2
    P_euler = (np.pi**2 * E * I) / (L**2)

    # Bounds for bisection
    P_min = 1e-6
    P_max = 0.99 * P_euler

    # Numerical solution
    P_critical = bisect(f, P_min, P_max, xtol=1e-6)

    return P_critical
