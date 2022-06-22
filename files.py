import os

from anyio import run_async_from_thread


files = os.listdir('.')

string = ""
for file in files[:4]:
    string += file + ' '
    print(file)

runnable_files = []
for file in files:
    if file != "files.py":
        runnable_files.append(file)

run_string = ""
for file in runnable_files[:50]:
    run_string += ' ' + file + ' '
print(run_string)


os.popen("ds9 -tile 56920a121-w2-int-1b_ra224.49725_dec-14.30182_asec600.000.fits  57032a115-w1-int-1b_ra229.527303_dec-12.54849_asec600.000.fits  56861b148-w2-int-1b_ra221.974927_dec-15.1014_asec600.000.fits  56972a118-w1-int-1b_ra226.807034_dec-13.52231_asec600.000.fits  56916a121-w2-int-1b_ra224.32138_dec-14.35936_asec600.000.fits  57004a118-w2-int-1b_ra228.250789_dec-13.01277_asec600.000.fits  57080a114-w1-int-1b_ra231.744228_dec-11.7121_asec600.000.fits  56997b142-w2-int-1b_ra227.978717_dec-13.11011_asec600.000.fits  57085b138-w2-int-1b_ra232.023706_dec-11.60406_asec600.000.fits  56901b146-w2-int-1b_ra223.708162_dec-14.55789_asec600.000.fits  56821b150-w2-int-1b_ra220.269578_dec-15.61039_asec600.000.fits  57065b139-w2-int-1b_ra231.094122_dec-11.9612_asec600.000.fits  56828a124-w2-int-1b_ra220.523523_dec-15.53626_asec600.000.fits  57060a114-w2-int-1b_ra230.816241_dec-12.06676_asec600.000.fits  56856a123-w2-int-1b_ra221.717232_dec-15.18_asec600.000.fits  56853b148-w2-int-1b_ra221.63165_dec-15.20592_asec600.000.fits  56953b144-w2-int-1b_ra226.002362_dec-13.79891_asec600.000.fits  57073b139-w2-int-1b_ra231.465281_dec-11.81937_asec600.000.fits  57117b137-w1-int-1b_ra233.523511_dec-11.0143_asec600.000.fits  56944a120-w1-int-1b_ra225.557631_dec-13.94949_asec600.000.fits  57081b139-w1-int-1b_ra231.837393_dec-11.67613_asec600.000.fits  56777b151-w2-int-1b_ra218.425981_dec-16.13127_asec600.000.fits  56990a054-w2-int-1b_ra227.617041_dec-13.23853_asec600.000.fits  57125b136-w2-int-1b_ra233.900486_dec-10.8636_asec600.000.fits  56785b151-w2-int-1b_ra218.758728_dec-16.03952_asec600.000.fits  57117b136-w2-int-1b_ra233.523318_dec-11.01444_asec600.000.fits  56988a117-w1-int-1b_ra227.526739_dec-13.27044_asec600.000.fits  56865b148-w1-int-1b_ra222.147072_dec-15.04857_asec600.000.fits  56877b147-w2-int-1b_ra222.664825_dec-14.88822_asec600.000.fits  56856a124-w1-int-1b_ra221.717408_dec-15.17989_asec600.000.fits  56845b148-w2-int-1b_ra221.289323_dec-15.30917_asec600.000.fits  56976a118-w1-int-1b_ra226.986476_dec-13.45993_asec600.000.fits  56824a125-w2-int-1b_ra220.354214_dec-15.58572_asec600.000.fits  57017b141-w2-int-1b_ra228.887496_dec-12.7828_asec600.000.fits  56956a120-w1-int-1b_ra226.091545_dec-13.7685_asec600.000.fits  56761b152-w2-int-1b_ra217.764227_dec-16.3106_asec600.000.fits  56865b147-w2-int-1b_ra222.146896_dec-15.04868_asec600.000.fits  56818a055-w1-int-1b_ra220.100518_dec-15.65945_asec600.000.fits  57133b135-w2-int-1b_ra234.278538_dec-10.71139_asec600.000.fits  56781b151-w1-int-1b_ra218.592213_dec-16.08556_asec600.000.fits  56773b152-w2-int-1b_ra218.260203_dec-16.17655_asec600.000.fits  57032a116-w1-int-1b_ra229.52749_dec-12.54836_asec600.000.fits  57076a113-w1-int-1b_ra231.558077_dec-11.78379_asec600.000.fits  56940a062-w2-int-1b_ra225.380316_dec-14.00904_asec600.000.fits  57013b141-w1-int-1b_ra228.705155_dec-12.84902_asec600.000.fits  56828a124-w1-int-1b_ra220.523523_dec-15.53626_asec600.000.fits  57060a114-w1-int-1b_ra230.816241_dec-12.06676_asec600.000.fits  57044a115-w2-int-1b_ra230.07828_dec-12.34411_asec600.000.fits  56953b144-w1-int-1b_ra226.002362_dec-13.79891_asec600.000.fits  57050a042-w1-int-1b_ra230.35451_dec-12.2408_asec600.000.fits")
