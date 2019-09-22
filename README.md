# quantum-annealing
d-waveの量子アニーリングマシンを使って遊ぶ

# 使い方
## アクティブなノードとエッジを確認する。
1. dwaveのtoken, solverを設定しておく。
~~~
> dwave config create
~~~
2.  アクティブなノード、エッジを標準出力する。
~~~
> python active.py
~~~
## ループを見つける。
~~~
> python find_loop.py
~~~
※モデルはあってるはずだが、最適解を見つけれない。
　マシンの問題？
