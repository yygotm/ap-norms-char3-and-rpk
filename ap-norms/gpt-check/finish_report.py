"""Validate our completed logs and atomically publish the requested report."""
import hashlib
import json
import re
from datetime import datetime, timedelta, timezone
from pathlib import Path

work = Path(__file__).resolve().parent
repo = work.parent.parent
pairs = [(3, 6), (3, 7), (3, 8), (3, 10), (3, 12), (9, 6), (9, 8), (27, 6)]
expected = [67, 157, 423, 2979, 22217, 355472, 21530128, 1549845936]

def objects(path):
    text = path.read_text(encoding='utf-8')
    assert not re.search(r'^(?:FAIL|ERROR)\b', text, re.M), path
    return [json.loads(s) for s in text.splitlines() if s.startswith('{')]

rows = []
for pair, count in zip(pairs, expected):
    name = f'{pair[0]}_{pair[1]}'
    obj = objects(work / f'{name}-full.log')
    entries = {}
    for mode in ('header', 'full', 'coverage', 'complete'):
        found = [x for x in obj if x['mode'] == mode]
        assert len(found) == 1, (name, mode, len(found))
        entries[mode] = found[0]
    h, f, c, end = (entries[m] for m in ('header', 'full', 'coverage', 'complete'))
    assert all((x['q'], x['n']) == pair for x in (h, f, c, end))
    assert h['polynomial_primitive'] is True
    assert h['witnesses'] == f['checked'] == count
    assert h['records'] == f['records'] == c['records']
    assert f['failed'] == 0 and f['failure_codes'] == [0] * 6
    assert f['threads'] == 8
    assert all(c[k] == 0 for k in ('invalid_j', 'duplicate_orbits', 'uncovered_elements', 'uncovered_orbits'))
    rows.append(entries)

batch = (work / 'batch.log').read_text(encoding='utf-8')
assert 'UTC batch finished:' in batch
assert not re.search(r'^(?:FAIL|ERROR)\b', batch, re.M)
assert len([x for x in objects(work / 'batch.log') if x['mode'] == 'complete']) == 7
assert objects(work / 'fixtures.log') == [{'mode': 'fixtures', 'passed': 232}]
jst = timezone(timedelta(hours=9))

def stamp(label):
    value = re.search(label + r': (\S+)', batch).group(1)
    return datetime.fromisoformat(value.replace('Z', '+00:00')).astimezone(jst)

begin, end = stamp('UTC batch started'), stamp('UTC batch finished')
hashes = {}
for name in ('check.cpp', 'check', 'fixtures.txt'):
    digest = hashlib.sha256((work / name).read_bytes()).hexdigest()
    assert digest in batch, ('hash differs from run-start snapshot', name)
    hashes[name] = digest

sample_sets = {}
for tag, filename in [('v1', '27_6-sample-v1.log'), ('final', '27_6-sample.log')]:
    samples = [x for x in objects(work / filename) if x['mode'] == 'sample']
    assert {x['region'] for x in samples} == {'head', 'tail', 'random'}
    assert all(x['checked'] == 30000 and x['failed'] == 0 and x['threads'] == 8 for x in samples)
    sample_sets[tag] = samples

table = []
for pair, e in zip(pairs, rows):
    h, f, c, done = (e[m] for m in ('header', 'full', 'coverage', 'complete'))
    table.append(f"| ({pair[0]},{pair[1]}) | {h['records']:,} | {h['witnesses']:,} | {f['checked']:,} | 0 | 0 | {f['seconds']:.6f} | {c['seconds']:.6f} | {done['seconds']:.6f} | PASS |")

sample_table = []
for tag in ('v1', 'final'):
    for s in sample_sets[tag]:
        sample_table.append(f"| {tag} | {s['region']} | {s['checked']:,} | 0 | {s['load_seconds']:.6f} | {s['check_seconds']:.6f} | {s['projected_full_seconds']/3600:.3f} |")

factor_table = []
for pair, e in zip(pairs, rows):
    q, n = pair
    M, v, factors = q**n-1, q**n-1, []
    for p in e['header']['prime_divisors']:
        a = 0
        while v % p == 0:
            v //= p
            a += 1
        factors.append(str(p) if a == 1 else f'{p}^{a}')
    assert v == 1
    factor_table.append(f"| ({q},{n}) | {M:,} | {' × '.join(factors)} | PASS |")

now = datetime.now(jst).strftime('%Y/%m/%d %H:%M:%S JST')
text = f'''# 候補31：独立 C++ 検査器による全8組の全件検査

作成：{now}

## 結論

**8組すべて全件 PASS。合計 {sum(expected):,} 証人を検査し、不合格0、未被覆軌道0、軌道代表の重複0。抜き取りだけの判定ではない。**

先の Python 版報告では (3,6) 以外を全件未実行としていたが、本報告では全8組について C++/OpenMP 8スレッドで全件検査を完了した。先の報告の未実行欄は、本報告の結果で更新される。

根拠は各組の `header` / `full` / `coverage` / `complete` レコードとバッチの正常終了。集計プログラム `finish_report.py` で総数・検査数の一致、全エラー数0、全完了記録、実行開始時と終了時のソース・実行ファイルのハッシュ一致を確認した。

## 組別結果

| 組 | レコード数 | 証人総数 | 全件検査数 | 不合格 | 未被覆軌道 | 証人検査 秒 | 軌道被覆 秒 | 合計 秒 | 判定 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
{chr(10).join(table)}

各組の合計秒は C++ の `complete.seconds`。ヘッダー検証・配列初期化など、`full()` に入る前の時間は含まない。全8組とも、不正な j、軌道代表の重複、未被覆の指数の数も0。

(3,6) の最終版 C++ の照合は2026/10/9 21:05:51–21:05:52 JST。残り7組のバッチは **{begin.strftime('%Y/%m/%d %H:%M:%S')}–{end.strftime('%Y/%m/%d %H:%M:%S')} JST**、実時間 **{(end-begin).total_seconds():,.0f}秒（約{(end-begin).total_seconds()/60:.2f}分）**。

## 独立性・実装

仕様は依頼書から取得した。作者側の実装・ログ・メモは読んでいない。証明書の保存場所は `<certificate dir>/<q>_<n>\\`。証明書データは変更せず、コミットもしない。

検査器は `ap-norms/gpt-check/check.cpp`。WSL の g++ 15.2.0、`-O3 -march=native -fopenmp -std=c++20` でビルドし、OpenMP を最大8スレッドに固定した。8組を同時に実行せず、順番に検査した。

F₃ の係数1と2を2つのビット平面に格納し、多項式の加算、乗算、f による剰余、二進法のべき乗を自作した。これは先の Python 版と同じ数学的な演算方法であり、相手側の計算コード・検査器・ログは参照していない。`poly_high_to_low` は文字列を逆順にして x⁰ からの係数と解釈した。

各 M は試し割りで独立に因数分解し、x^M=1、および各素因数 p に対する x^(M/p)≠1 を確認した。3^N 個の元を持つ商環に位数3^N−1の単元が存在するので、非零元はすべてその冪であり、商環が体であることも従う。

| 組 | M | 素因数分解 | x の位数 M |
|---|---:|---|---|
{chr(10).join(factor_table)}

各証人について α=x^i、β=x^j、α+β を上記の演算で作り、両項の原始性を素因数ごとのべき乗条件で、両ノルムを M/(q−1) 乗で確認した。ノルムの期待値は g_q^c を自作演算で求め、原始指数 c を昇順に並べた行優先のスロットと照合した。

**正規性は F₃ 上の階数で判定した。** F_q の F₃ 基底を 1,g_q,…,g_q^(k−1) とし、これらと q-共役 α,α^q,…,α^(q^(n−1)) の全積から N=kn 本の F₃ ベクトルを作る。この N×N 行列が階数 N であることは、元の共役が F_q 上の基底であることと同値である。α または α+β の片方が正規なら条件を満たす。相手側の正規性判定方法は読んでいないので、「相手と違う方法だった」とは主張しない。

高速化は、同じ元のべき乗で中間の二乗列を再利用し、α=x^i の確定済みの原始性・ノルム・正規性を指数 i ごとに保存するもの。キャッシュは組ごとに新規作成し、初期値は未計算。未計算の元は上記の独立演算で検査し、値を atomic な1バイトに保存する。再利用しても β や指定ノルムのスロットの検査を省略しない。α+β は各証人で構成して検査する。

自作 Python 版との照合では、8組合計232例（0、1、2、x、最大の密な整数表現、および固定seedによる無作為元）について、加算・乗算・べき乗・正規性が一致した。共有二乗列によるべき乗も通常のべき乗と一致した。(3,6) は C++ で67件・不合格0・未被覆軌道0となり、Python の全件結果と一致した。

軌道被覆は各 j から 3^t j と 3^t j+M/2（mod M、0≤t<N）を列挙してビット集合に記録し、全 M 個の指数について被覆を確認した。未被覆があればさらにその軌道を列挙して未被覆軌道数を数える。全8組とも0だった。

## 速度測定・実行の承認

(27,6) は `cert_*.bin` 16本、合計 **6,242,435,020バイト**であることを確認してから着手した。先頭30,000証人、末尾30,000証人、それらを除く範囲から無作為30,000証人（重複なし、seed=202610092056）を測定した。v1 と最終版で同じ選択方法を使い、どちらも90,000証人の不合格0だった。

| 版 | 標本 | 検査数 | 不合格 | 読み込み 秒 | 検査 秒 | (27,6)全件外挿 時間 |
|---|---|---:|---:|---:|---:|---:|
{chr(10).join(sample_table)}

v1 の8スレッド測定では (27,6) 単独で約3.97–4.39時間と見積もられた。最終版では同じ数学的判定の計算結果を再利用して短縮した。上表は短い標本の速度を外挿した値であり、全件実行時間とは区別する。速度差の原因を分解して測定したわけではない。

## 数学側の判断

前報告の判断を維持する。β→β³ は原始性と F_q 上の正規性を保ち、ノルムの組を c→c³ で全単射に置換する。β→−β は2項の交換で対応する。したがって全順序付き原始ノルム対を検査する本証明書の軌道還元は正しい。

q=3 では F₃ の原始元は −1 のみであり、拡大体の原始元のノルムは必ず −1。よって指定ノルム条件は自動的に満たされる。[Lemos–Neumann–Ribas, Theorem 1.1(i), Table 6](https://arxiv.org/pdf/2208.02876) は任意の非零 β に対する原始・片方正規な2項の存在を、表6の候補以外で保証し、(3,9) と (3,16) は表6にない。したがってこの2組も 𝒮₂ に入る。これらはノルム指定なしの定理に q=3 のノルムの観察を組み合わせた判断である。原典確認の詳細は前報告を参照。

[Chatterjee–Sharma–Shukla–Tiwari, Theorem 5.7](https://arxiv.org/pdf/2401.01819) の残る候補10組について、8組は本全件検査、2組は上記の既存定理により条件を満たすという結論になる。ただし、この報告で論文全体の証明を再査読したわけではない。

## 再現用ファイルとハッシュ

作業場所：`ap-norms/gpt-check/`

- `check.cpp`：最終 C++ 検査器。
- `check-v1.cpp`：最初の C++ 速度測定時の実装。
- `verify.py`：独立に作成した Python 版。
- `make_fixtures.py`, `fixtures.txt`, `fixtures.log`：Python/C++ の232例照合。
- `run.sh`：ビルド・単独実行・残り7組のバッチ実行。
- `*-full.log`, `batch.log`：全件の検査記録。
- `27_6-sample-v1.log`, `27_6-sample.log`：速度測定の記録。
- `finish_report.py`：完了記録・件数・ハッシュを照合して本報告を生成。

```text
SHA256 check.cpp    {hashes['check.cpp']}
SHA256 check        {hashes['check']}
SHA256 fixtures.txt {hashes['fixtures.txt']}
```

単独再実行例：

```powershell
DATA=<certificate dir> bash ap-norms/gpt-check/run.sh full 27_6 8
```

再現には上記データ保存場所の証明書と WSL の C++/OpenMP 環境が必要。コードはローカルに置き、コミット・外部送信はしていない。
'''

output = repo / 'reports/20261009-candidate31-independent-check-part2.md.regenerated'
temp = output.with_name(output.name + '.tmp')
temp.write_text(text, encoding='utf-8', newline='\n')
temp.replace(output)
print(json.dumps({'output': str(output), 'witnesses': sum(expected),
                  'all_pass': True, 'batch_wall_seconds': (end-begin).total_seconds()}, ensure_ascii=False))
