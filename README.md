# PPGR - 1. seminarski
Ovaj repozitorijum je namenjen potrebama kursa Primena projektivne geometrije u racunarstvu na Matematickom fakultetu.

Zadatak: 
Računanje projektivnog preslikavanja i otklanjanje distorzije.
Implementirani:
1) Naivni algoritam
2) DLT algoritam (i poredjenje)
3) Modifikovani DLT algoritam (sa normalizacijom)
4) Aplikacija koja ucitava bmp sliku, korisnik bira (mišem ili unosi koordinate) 4 piksela koji se slikaju u pravougaonik, a aplikacija vraca sliku sa otklonjenom projektivnom distorzijom

Jezik: C++

Biblioteke: Eigen, CImg

Kod se na Linux-u iz terminala prevodi sa: g++ -std=c++17 -o out code.cpp -lm -lpthread -lX11, a pokrece sa ./out
