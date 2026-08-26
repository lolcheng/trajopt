from core.variable import Variables, VarBlock
import numpy as np

if __name__ == "__main__":
    n_ctrl = 5
    vx = VarBlock(name="Px", shape=(n_ctrl,), lb=-10, ub=10, scale=1.0, init=np.linspace(0, 1, n_ctrl))
    vy = VarBlock(name="Py", shape=(n_ctrl,), lb=-5, ub=5, scale=2.0, init=0.0)
    vs = Variables([vx, vy])
    print(vs.summary())
    x0 = vs.initial_x()
    print("x0 (scaled):", x0)
    parts = vs.unpack(x0)
    print("Px unscaled:", parts["Px"]) 
    print("Py unscaled:", parts["Py"]) 
    assert vs.is_feasible(x0)