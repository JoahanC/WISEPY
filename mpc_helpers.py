_mpc_hex = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"

def unpackMPCname(compact):

    if len(compact) == 5:

        if compact[4] == "P":

            # Periodic comet

            outn = "{:s}P".format(str(int(compact[0:4])))

        elif compact[4] == "S":

            # Natural Satellite

            outn = compact

        elif compact[0] == "~":

            # Numbered object at or above 620000

            dig1 = _mpc_hex.index(compact[1]) * (62**3)

            dig2 = _mpc_hex.index(compact[2]) * (62**2)

            dig3 = _mpc_hex.index(compact[3]) * 62

            dig4 = _mpc_hex.index(compact[4])

            num = 620000 + dig1 + dig2 + dig3 + dig4

            outn = "({:s})".format(str(num))

        else:

            # Numbered object

            outn = "({:s})".format(

                str(_mpc_hex.index(compact[0]) * 10000 + int(compact[1:]))

            )

    elif len(compact) == 7:

        # PLS object

        if compact[0:3] in ["PLS", "T1S", "T2S", "T3S"]:

            outn = "{:s} {:s}-{:s}".format(compact[3:], compact[0], compact[1])

        else:

            # Unnumbered asteroid

            year = int(_mpc_hex.index(compact[0]) * 100) + int(compact[1:3])

            if compact[4:6] == "00":

                prov = "{:1s}{:1s}".format(compact[3], compact[6])

            else:

                prov = "{:1s}{:1s}{:d}".format(

                    compact[3],

                    compact[6],

                    int(_mpc_hex.index(compact[4])) * 10 + int(compact[5]),

                )

            outn = "{:s} {:s}".format(str(year), prov)

    elif len(compact) == 8:

        # parabolic comet

        year = _mpc_hex.index(compact[1]) * 100 + int(compact[2:4])

        if compact[7] not in _mpc_hex[0:10]:

            prov = "{:s}{:s}{:s}".format(

                compact[4],

                compact[7],

                str(_mpc_hex.index(compact[5]) * 10 + int(compact[6])),

            )

            outn = "{:s}/{:s} {:s}".format(compact[0], str(year), prov)

        else:

            outn = "{:s}/{:s} {:s}{:s}".format(

                compact[0], str(year), compact[4], str(int(compact[5:7]))

            )

    else:

        raise NotImplementedError(

            "This designation could not be unpacked {:s}".format(compact)

        )

    return outn

 

 

def conv_num2mpcformat(n):

    """

    TBD

    :param n:

    :return:

    """

    nn2c = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"

    if n < 100000:

        return "%05d" % n

    elif n < 620000:

        nn = int(n / 10000.0)

        c = nn2c[nn]

        return "{:1s}{:04d}".format(c, n - nn * 10000)

    else:

        # For numbers larger than 620,000 the MPC has defined a new packing

        # scheme. Code by J. Masiero

        nn = n - 620000

        dig4 = int(nn % 62)

        hold3 = (nn - dig4) / 62.0

        dig3 = int(hold3 % 62)

        hold2 = (hold3 - dig3) / 62.0

        dig2 = int(hold2 % 62)

        hold1 = (hold2 - dig2) / 62.0

        dig1 = int(hold1 % 62)

        return "~{:1s}{:1s}{:1s}{:1s}".format(

            nn2c[dig1], nn2c[dig2], nn2c[dig3], nn2c[dig4]

        )

 

 

def conv_des2compact(y, des):

    """

    TBD

    :param n:

    :return:

    """

    if des == "P-L":

        return "PLS{:04d}".format(int(y))

    elif des == "T-3":

        return "T3S{:04d}".format(int(y))

    elif des == "T-2":

        return "T2S{:04d}".format(int(y))

    elif des == "T-1":

        return "T1S{:04d}".format(int(y))

    nn2c = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"

    y = int(y)

 

    # if y is not in the range of [1800,2100), an error would occur

    if 1800 <= y < 1900:

        yy = "I{:02d}".format(y - 1800)

    elif 1900 <= y < 2000:

        yy = "J{:02d}".format(y - 1900)

    elif 2000 <= y < 2100:

        yy = "K{:02d}".format(y - 2000)

 

    fl = des[0]

    sl = des[1]

    if len(des) > 2:

        num = int(des[2:])

        if num > 99:

            i = int(num / 10.0)

            mid = "{:1s}{:1d}".format(nn2c[i], num - i * 10)

        else:

            mid = "{:02d}".format(num)

        return "{:3s}{:1s}{:2s}{:1s}".format(yy, fl, mid, sl)

    else:

        return "{:3s}{:1s}00{:1s}".format(yy, fl, sl)