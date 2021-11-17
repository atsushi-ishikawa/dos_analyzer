## VASPでDOSを自動的に計算してmachine-learning

* 流れ
1. VASPの計算を行う: **adsorption.py**

* bimetallic alloyを想定。上記pythonスクリプトの引数で計算対象を指定する。例えば、
`python adsorption.py Pt` --> pure metalを計算
`python adsorption.py Pt Pd 50` --> Pt0.5Pd0.5のalloyを計算。第三引数は第一引数の元素のパーセンテージ。
* adsorption.pyを実行すると計算条件と計算結果がjsonファイル(**surf_data.json**)に追記される。汎関数や吸着エネルギーなど。
* *DOSCAR* と *vasprun_{system}.xml* が出力される。systemはalloy surfaceの名前(Ru0.5Rh0.5_111など)
* bimetallic alloyの成分を変えながら、一連の計算をsubmitするには**submit_ads.py**を使う
* submit_ads.py内でadsorption.pyにジョブを投げるようになっている。したがって計算機センター等では`python submit_ads.py`を実行すると一連の計算が始まる

2. DOSCARのpeakを検出する: **analyze_dos.py**
* DOSCARを読み込んでs/p/d-DOSのピークを検出
* 計算データの入ったjsonファイル("surf_data.json")が必要
* 上記jsonファイルとは別のjsonファイルに結果が出力される(*tmpout.json*)。ここにDOSのピークの情報と吸着エネルギーが書かれる <-- 改善の余地あり
* １つのDOSCARを解析する場合
```bash
python analize_dos.py [system_name] [orbital(s/p/d)]
 example: python analize_dos.py Rh0.2Cu0.8_111 d
```
* フォルダ内の全てのDOSCARに行う場合
```bash
python all_analize.py
```

3. 検出したピークを記述子にして吸着エネルギーの回帰を行う: **regression.py**
* analyze_dos.pyから出力されたtmpout.jsonを使う
* sklearnを利用。現状ではOLS (ordinary linear regression), Ridge, LASSO, random forestを行う
* BOLASSO (bootstrapped LASSSO)はRにしかライブラリがないので**do_bolasso.r**で実行する。regression.pyからはこのRのスクリプトを呼ぶようになっている
* あああ

---

### その他スクリプト
* **plotdos.py** : DOSCARをプロット
```bash
python plotdos.py Rh0.2Cu0.8_111 ; open DOS_Rh0.2Cu0.8_111.png
```
* **waterfall_plot.py** : DOSCARが縦に一斉に並んだようなplotを作る
* **colormap.py** : alloyの要素に対して吸着エネルギーのカラーマップを作る
