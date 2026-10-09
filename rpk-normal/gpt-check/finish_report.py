from pathlib import Path
import hashlib,json,math,re,datetime
ROOT=Path(__file__).resolve().parent
import os
DATA=Path(os.environ.get('DATA', str(ROOT.parent / 'witnesses')))
DEST=ROOT.parents[1]/'reports'/'20261010-candidate32-independent-check.md.regenerated'  # does not overwrite the archived report
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def parse(name):
    text=(ROOT/name).read_text();rows=[]
    for line in text.splitlines():
        if line.startswith('RESULT '):
            z=line.split();d=dict(x.split('=',1) for x in z[3:]);d.update(group=z[1],pair=z[2]);rows.append(d)
    return text,rows
def fmt(x):return f'{int(x):,}'
def run():
    full,rows=parse('all-full.log');bench,br=parse('all-bench.log')
    assert 'COMPLETE mode=full pair=all' in full and 'UTC completed:' in full and 'ERROR' not in full
    assert len(rows)==78 and len(br)==78
    aa=[r for r in rows if r['group']=='A'];bb=[r for r in rows if r['group']=='B']
    assert len(aa)==21 and len(bb)==57 and len({r['pair'] for r in rows})==78
    for r in aa:
        assert r['mode']=='full' and r['scanned']==r['potential'] and int(r['bad'])==0
        assert int(r['missing'])==(1 if r['pair']=='2_8' else 0)
        m=int(r['pair'].split('_')[1]);assert int(r['orbits'])*m==int(r['primitive'])
        line=next(x for x in full.splitlines() if x.startswith('COUNTS A '+r['pair']+' '))
        fields=dict(x.split('=') for x in line.split()[3:]);counts=[int(x) for x in fields['cells'].strip(',').split(',')]
        traces=[int(x) for x in fields['traces'].strip(',').split(',')]
        assert len(counts)==int(r['cells']) and sum(counts)==sum(traces)==int(r['normal'])
        assert sum(x==0 for x in counts)==int(r['missing'])
    for r in bb:
        assert r['mode']=='full' and r['checked']==r['total']
        assert all(int(r[k])==0 for k in ('bad','order_bad','rank_bad','cell_bad','missing','duplicates'))
    manifest=json.loads((ROOT/'input-manifest.json').read_text())
    for x in manifest:
        assert sha(DATA/(x['pair']+'.hdr'))==x['header_sha256']
        assert sha(DATA/(x['pair']+'.bin'))==x['binary_sha256']
    assert sum(int(r['checked']) for r in bb)==462321
    ref=json.loads((ROOT/'reference-results.json').read_text())
    assert ref['B_checked']==171 and ref['B_failed']==0 and ref['A_2_8_inverse_trace']=={'0':0,'1':8}
    for name in ('check.cpp','check','cases.txt'):assert sha(ROOT/name) in full
    estimate=sum(float(r['estimate_s']) for r in br);assert estimate<7200
    now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9)))
    durations={g:sum(float(r['total_s']) for r in rows if r['group']==g) for g in ('A','B')}
    wall=float(re.search(r'COMPLETE mode=full pair=all elapsed_s=([\d.]+)',full).group(1))
    sums={k:sum(int(r[k]) for r in aa) for k in ('field','potential','primitive','orbits','normal','cells','missing')}
    out=[]
    def put(s=''):out.append(s)
    put('# 候補32 独立検査結果')
    put(f'\n作成: {now:%Y-%m-%d %H:%M:%S} JST。依頼: 2026-10-10 04:05 の candidate32-independent-checker。')
    put('対象: r=3, k=1 の指定78組。検査器は独自の C++20 / OpenMP、最大8スレッド。')
    put('\n## 結論')
    put('- **A: 21組の全数分類が完了。20組で許される全セルを被覆し、(2,8) だけ (a,b)=(0,1) が空であることを確認。依頼書の計算上の主張と一致。**')
    put('- **B: 57組・462,321証人を全件検査し、全件合格。既約性・位数・1-normal性・係数セル・被覆の不一致は0、重複0、未被覆0。**')
    put('- ノルムの正しい位数条件は ord(b)=(q−1)/gcd(3,q−1)。Lemma 3.4 の w-free という記述はこの条件と一致せず、Theorem 3.7 の任意の非零 b という記述はそのままでは偽。これは w=3 だけでなく w=1 の場合にも制限を与える。')
    put('- **この検査が保証する有限個の存在結果は、b を上記の位数に制限した問題に対するもの。論文の無制限な定理を検証・修復したという意味ではない。**')
    put('\n## 数学側の独立確認')
    put('原文は [arXiv:2501.04999 v3, 2025-04-16](https://arxiv.org/pdf/2501.04999) の定義（p.2）、Lemma 3.1・3.4（p.5）、Theorem 3.7 とその直前（p.6）を確認した。以下はその記述に対する独立の群論的検査である。')
    put('\nQ₀=q^m−1, d=q−1 と書き、ε を F_{q^m}^* の生成元、ξ=ε^t とする。ord(ξ)=Q₀/3 と gcd(t,Q₀)=3 は同値である。N(ε)=ε^{Q₀/d} は位数 d を持ち、d|Q₀ より')
    put('\n```text\ngcd(t,d) = gcd(gcd(t,Q₀),d) = gcd(3,d) = w,\nord(N(ξ)) = d/gcd(t,d) = (q−1)/w.\n```')
    put('\nしたがって依頼書の位数の導出は正しい。逆に、位数 d/w の元は原始元の3乗として表せる（巡回群の指数を素因数ごとに持ち上げる）。この条件は「3-primitive 元のノルムとして可能」の必要十分条件であり、1-normal 条件は別途検査が必要である。')
    put('\n原文の l-free の定義は gcd(l,(q−1)/ord(b))=1。q=4 では w=3 であり、可能なノルムは b=1 だけだが、これは3-freeではない。一方、F₄^* の位数3の元は3-freeだが、3-primitive 元のノルムにはなれない。Lemma 3.4 はこの例で両方向とも成立しない。q=5 では w=1 なので全非零元が1-freeだが、実際に可能なノルムは位数4の元だけであり、b=1 などは不可能である。')
    put('\nTheorem 3.7 の直前には、原始元 c を c^r=b となるように選ぶという段階がある。この選択が可能な b には上記の制約があるため、定理文の任意の非零 b へ広げることはできない。')
    put('\n### 定理の印刷された十分条件も満たす具体的反例')
    put('q=4, m=64, r=3, k=1, g=X−1 とする。g|X^64−1 かつ k<m/2。b を F₄^* の位数3の元とし、例えば a=0 とする。3-primitive 元のノルムは必ず1なので、求める ξ は存在しない。')
    put('\n一方、独立に積と各因子の素数性を確認した分解は')
    put('\n```text\n4^64−1 = 3 · 5 · 17 · 257 · 641 · 65537 · 274177 · 6700417 · 67280421310721.\n```')
    put('\n原文の Q（ここでの Q₀ とは別記号）は Q₀ のうち q−1=3 と互いに素な最大の約数だから W(Q)=2^8。標数2では (X^64−1)/(X−1)=(X−1)^63 なので、その W は2。よって印刷された不等式は')
    put('\n```text\nq^(m/2−k−2) = 4^29 = 288230376151711744 > 3 · 256 · 2 = 1536\n```')
    put('\nと厳密に成立する。従って、これは不等式の適用範囲外の例を挙げただけではなく、Theorem 3.7 の文面に対する反例である。g を割り切る条件の省略など、論文の他の論点の全面査読は本依頼の範囲外。')
    put('\n### 末尾係数と最小多項式の次数')
    put('共役を y₀,…,y_{m−1} とすると、b=∏yᵢ、a=Σᵢ∏ⱼ≠ᵢyⱼ=b·Tr(ξ⁻¹) である。C++ は積と第(m−1)基本対称式を漸化式で計算し、Python 照合では逆元・トレース・ノルムを直接計算した。m≥7 で1-normalなら共役の張る空間は m−1 次元であり、真部分体の次数は m/2 以下なので ξ の最小多項式次数は m。係数の符号は依頼書の a_{m−1}, a_m の規約で統一した。')
    put('\n## 検査方法と完全性')
    put('作者側の実装・ログ・メモは読んでいない。参照した入力は指定依頼書、論文PDF、提供された .hdr/.bin と witnesses.done。実装・ログ・補助資料は rpk-normal/gpt-check/ のみ。証人は変更せず、実行前後のSHA256の一致を57組すべてで確認した。コミットしていない。')
    put('\nQ₀ の素因数分解は prepare.py が Pollard rho と整数演算で独立に求めた。各素因数は2^64未満で、決定的 Miller–Rabin の7基底 {2,325,9375,28178,450775,9780504,1795265022} を使用。C++ 側でも素数性と積=Q₀ を再検査した。Q₀ と冪指数には unsigned __int128 を使用し、(883,7) などの64ビット超過を避けた。分解結果は cases.txt に保存。')
    put('\n多項式既約性は Rabin 判定：X^{p^N}=X mod f、および各素数 ℓ|N に対して gcd(f,X^{p^{N/ℓ}}−X)=1。Fp 上の加減乗算・多項式剰余・gcd を実装した。2進体はビット演算、奇標数では係数配列を使用。')
    put('\n### A: 全体の走査と共役による重複計算の省略')
    put('各組について独自の既約多項式 f を選び、X の位数が Q₀ であることを全素因数で検査した。1,X,…,X^{Q₀−1} をすべて生成するので、非零元全体を一度ずつ走査する。0は3-primitiveではない。位数条件から指数が3の倍数でない元は厳密に除外し、ξ=X^{3t}, 0≤t<Q₀/3 を全走査、gcd(t,Q₀/3)=1 の元をすべて分類した。途中で被覆が揃っても打ち切っていない。')
    put('\n1-normal性とセルは q-Frobenius 共役で不変なので、同じ共役軌道では最小指数の元について一度だけ階数とセルを計算し、軌道長を件数に加算した。すべての該当軌道の長さが m であること、および「軌道数×m=φ(Q₀/3)=位数条件を満たす元の総数」を実行時に照合した。これは無作為標本検査ではなく、全元を覆う厳密な共役分類である。')
    put('\nFq の Fp 基底 u₀,…,u_{e−1} を取り、uⱼξ^{q^i}（0≤i<m, 0≤j<e）の Fp 上の階数が e(m−1) であることを判定した。階数がこの値を超えた場合だけ早期に非該当と判定。許される b と全 a を列挙してセル別件数を数えた。')
    put(f'\n総対象元数（0を含む）={fmt(sums["field"])}、生成した非零元={fmt(sums["field"]-21)}、指数 t の全検査数={fmt(sums["potential"])}、3-primitive 元={fmt(sums["primitive"])}、その共役軌道={fmt(sums["orbits"])}、3-primitive 1-normal 元={fmt(sums["normal"])}。期待265セル中264セルを被覆し、未被覆1セルは (2,8) の指定例外と一致。')
    put('\n### B: 証人とセルの全件検査')
    put('57組の modulus を既約と確認。各バイナリの長さ=2N×cells、係数の範囲、全レコード数を検査。各 ξ について ξ^{Q₀/3}=1 かつ各素数 ℓ|(Q₀/3) について ξ^{Q₀/(3ℓ)}≠1 を確かめた。Fq を q乗で固定される部分体として構成し、Fq 基底の Fp 独立性も検査。その基底を用いて階数 e(m−1) を検査した。')
    put('\na,b を共役から再計算し、それらが Fq に属すること、ord(b)=(q−1)/w、全セルの相異性、レコード数=q·φ((q−1)/w)、全期待セルの被覆を確認した。ヘッダーの cells を正しいと仮定せず、独立計算した期待数と照合した。')
    put('\n### 別方式の照合')
    put('57組それぞれ先頭・中央・末尾の3件、合計171件を別の Python 実装で再計算した。k-normal性は拡大体上の gcd(X^m−1, ξX^{m−1}+ξ^qX^{m−2}+…+ξ^{q^{m−1}}) の次数で判定し、C++ の行列階数とは別の判定法を使った。位数、係数 a,b、k=1 の全171件が一致。さらに (2,8) は256個の係数表示を直接全列挙し、冪表・軌道短縮を用いずに同じ分布を確認した。')
    put('\n## (2,8) の完全な分布')
    put('\n| Tr(ξ⁻¹)=a | b | 3-primitive 1-normal 元の数 |\n|---:|---:|---:|\n| 0 | 1 | 0 |\n| 1 | 1 | 8 |')
    put('\n使用体は F₂[X]/(X⁸+X⁷+X⁵+X³+1)。低次係数をビットにした整数表示で、該当8元は 15,44,63,85,92,109,163,219。位数はすべて85、F₂階数は7。これらは一つの長さ8の Frobenius 軌道をなす。3-primitive 元全体は64個。')
    put('\n## 組別結果')
    put('\nA表の「位数検査」は t の検査数であり、それ以前に冪表で全非零元を生成した。「3-primitive全数」はすべての軌道を展開した件数。「不整合」は検査器の前提・件数・セル判定の不一致数で、性質を満たさない一般の体元の数ではない。')
    put('\n| (q,m) | 全体元数 | 位数検査 | 3-primitive全数 | 1-normal該当数 | 被覆/期待セル | 不整合 | 秒 | 結果 |')
    put('|---|---:|---:|---:|---:|---:|---:|---:|---|')
    for r in aa:
        q,m=r['pair'].split('_');covered=int(r['cells'])-int(r['missing']);verdict='例外確認' if r['pair']=='2_8' else '合格'
        put(f'| ({q},{m}) | {fmt(r["field"])} | {fmt(r["scanned"])} | {fmt(r["primitive"])} | {fmt(r["normal"])} | {covered}/{r["cells"]} | {r["bad"]} | {float(r["total_s"]):.6f} | {verdict} |')
    put('\nBは全組で既約性検査に合格。以下の証人数・検査数・期待セル数はすべて等しく、不合格・重複・未被覆が0。')
    put('\n| (q,m) | 証人数=検査数=期待セル数 | 不合格 | 重複 | 未被覆 | 秒 | 結果 |')
    put('|---|---:|---:|---:|---:|---:|---|')
    for r in bb:
        q,m=r['pair'].split('_');put(f'| ({q},{m}) | {fmt(r["checked"])} | 0 | 0 | 0 | {float(r["total_s"]):.6f} | 合格 |')
    put('\n## 速度測定・実行記録')
    put(f'全78組の測定からの合計見積もりは {estimate:.6f} 秒（A {sum(float(r["estimate_s"]) for r in br if r["group"]=="A"):.6f} 秒、B {sum(float(r["estimate_s"]) for r in br if r["group"]=="B"):.6f} 秒）で、2時間の停止基準を下回ったため全件検査を実施した。Aは最大300,000指数（先頭・末尾・固定種の擬似乱数各100,000）、Bは組ごとに最大3,000件を均等間隔で測定。小さい組は全件。構築時間を含めて組ごとに外挿した。')
    put(f'\n本検査は **2026-10-10 04:29:56 JST〜04:30:17 JST**、全体 {wall:.6f} 秒で正常終了（終了コード0）。組別時間の合計は A {durations["A"]:.6f} 秒、B {durations["B"]:.6f} 秒。総実時間との差にはメモリ解放、照合用出力などが含まれる。Python照合は {ref["elapsed_s"]:.6f} 秒。これらは機械実行時間であり、実装・原文確認・報告作成の時間は含まない。')
    put('\n概算の「3万行程度」に対して、実ファイルは462,321レコード、バイナリ総量9,077,976バイトだった。件数の概算のみが異なり、証人内容の不一致はなかった。提供された全57組と指定78組からAを除いた集合の一致も確認した。')
    put('\nビルド: g++ (Ubuntu 15.2.0-16ubuntu1) 15.2.0、-O3 -march=native -fopenmp -std=c++20 -Wall -Wextra。omp_set_dynamic(0), omp_set_num_threads(8)。')
    put('\n```text\nDATA=rpk-normal/witnesses bash rpk-normal/gpt-check/run.sh bench all\nDATA=rpk-normal/witnesses bash rpk-normal/gpt-check/run.sh full all\n```')
    put('\n## 証拠ファイルと SHA256')
    put('基準ディレクトリ: `rpk-normal/gpt-check/`。all-full.log に全78組の結果・Aの各セル件数・独自 modulus、cases.txt に全分解、input-manifest.json に全入力のハッシュ、reference-results.json に別方式照合の結果を記録した。')
    put('\n| ファイル | SHA256 |\n|---|---|')
    names=['check.cpp','check','run.sh','prepare.py','cases.txt','input-manifest.json','all-bench.log','all-full.log','reference.py','reference.log','reference-results.json','finish_report.py']
    for name in names:put(f'| {name} | `{sha(ROOT/name)}` |')
    put('\n速度測定後、本検査の前に照合用3件の出力を追加したため測定時と本検査時のソース・実行ファイルのSHA256は異なる。上表は本検査に使用した版で、そのハッシュが all-full.log 冒頭と一致することを確認済み。判定・列挙の計算部分は同一。')
    put('\nこれは独立実装による再現可能な計算検査であり、形式証明支援系で検証した形式証明ではない。全78組を超えた範囲や、論文の篩証明全体を検証したとは主張しない。')
    text='\n'.join(out)+'\n';tmp=DEST.with_suffix('.md.tmp');assert not DEST.exists()
    tmp.write_text(text,encoding='utf-8');assert tmp.read_text(encoding='utf-8')==text;tmp.replace(DEST)
    print('REPORT',DEST);print('SHA256',sha(DEST));print('A=21 claims confirmed B=57 all pass records=462321 source_and_input_hashes_unchanged=true')
if __name__=='__main__':run()
