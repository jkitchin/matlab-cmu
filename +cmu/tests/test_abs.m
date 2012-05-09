function test_abs

u = cmu.units;

a = -1*u.kg;

assertEqual(abs(a),1*u.kg)

assertEqual(abs([-1 -1]*u.s),[1 1]*u.s)