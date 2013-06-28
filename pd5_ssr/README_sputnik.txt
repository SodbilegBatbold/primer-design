Our thanks to Chris Abajian for allowing the use of Sputnik in primrose_ssr.
http://espressosoftware.com/sputnik/

We have modified his code by converting sputnik.c to sputnik.cpp, adding a constructor and allowing the parameters (previously set by #defines) to be set in the constructor.

Then, towards the end of Sputnik's findRepeats function we create a 'primer_factory', an object of new class 'sputnik_primrose' that we created for this pd5_ssr application. We then call primer_factory.find_primers(...) to find primers suitable for the repeat that Sputnik has just discovered.

