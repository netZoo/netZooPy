from netZooPy import dragon
from netZooPy import smaug


def test_smaug():
    print('Start SMAUG run ...')
    n = 1000
    p1 = 500
    p2 = 100
    expr, meth, theta, _ = dragon.simulate_dragon_data(eta11=0.005, eta12=0.005, eta22=0.05,
                                                       p1=p1, p2=p2, epsilon=[0.1, 0.1], n=n, seed=123)
    smaug_obj = smaug.Smaug(expr, meth)
    smaug_obj.run_smaug(keep_in_memory=False, output_fmt='.txt')
