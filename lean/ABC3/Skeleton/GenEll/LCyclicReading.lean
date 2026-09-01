/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.Lemma37Full
import ABC3.Found.GenEll.EllModuliGalois
import ABC3.Meta.Claim

/-!
# 第 1152 ブロック —— **`l`-巡回部分群の 2 つの読み方が食い違っている**（`Skeleton`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.18–p.20。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

## ★★★★★★★★★★★★★★★★これは何か——**逸脱が消費側に響いた**

原文の `H_L ⊆ E_L` は **`l`-巡回部分群スキーム**である。本企画はこれを
**2 通りに読んで**しまっており、その 2 つが繋がっていない。

| 読み | 定義 | どこで使うか |
|---|---|---|
| (V) `HasLCyclicVelu` | ☆**`L` 有理点** `Q` で `addOrderOf Q = l`、その Vélu の商が楕円曲線で半安定 | `Lemma 3.5` → `Lemma 3.7`（第 1151） |
| (J) `HasLCyclicJ` | ☆`E[l]` の中の **`Gal`-安定な直線** | `Theorem 3.8`（`imageContainsSL2J_of_alpha'` の `hno`） |

★**含意は片側だけ立つ**: 有理点 `Q` は `⟨Q⟩` を点ごとに固定するので

    (V) ⟹ (J),   したがって   ¬(J) ⟹ ¬(V)

である。☆(J) ⟹ (V) は**偽**——安定な直線の点は `L` 有理とは限らない。

## ★★★★★★★★★★★★食い違いの現れ方

`Lemma 3.7` の第 3 の主張は「`l`-巡回なら `[E] ∈ Exc`」、
`Theorem 3.8` はその対偶「`[E] ∉ Exc` なら `l`-巡回でない」を使う。

| | `Lemma 3.7` が与えるもの | `Theorem 3.8` が要るもの |
|---|---|---|
| 第 1151 の形 | `¬ HasLCyclicVelu` | `¬ HasLCyclicJ` |

★★**向きが合わない**——`¬(J) ⟹ ¬(V)` であって逆ではない。
☆すなわち第 1151 の `lemma_3_7` は `Theorem 3.8` を**まだ養えない**。

★★★これは `CLAUDE.md` の「逸脱」の規約が名指ししている場合である
——『後続の証明に影響が出ないならば読み替えを許容する』。**影響が出た**ので記録する。

## ★★★★★★★★★★★★★★★★節点（進捗枠 **2.7 / 3**）

`Lemma 3.7` の第 3 の主張を **(J) の側で**言い直すのに要るもの:

| # | 節点 | 内容 | 重み |
|---|---|---|---|
| 1 | `veluQuotientFull_descends` | ★★**閉じた**（第 1153 の `fixesCoeffs_veluQuotientFull` ＋ 第 1154 の `mem_range_of_fixed`） | 12 → **0** |
| 2 | `lemma_3_5_height_ineq_stableLine` | ☆第 1 で作った `L` 上の商に対して `Lemma 3.5` の高さ不等式を回す。★**良い素点の側は第 1155-1158 で閉じた**（下記） | 10 → **1** |
| 3 | `lemma_3_7_stableLine` | ★★**第 1166 で取れた**——`(†)` を仮説 `hdag` に切り出した形で。☆**残るギャップは `hdag` ただ 1 つ** | 4 → **0** |

☆総重み 26 → **4**（第 1153-1154 で節点 1、第 1155-1158 で節点 2 の良い素点の側、
第 1166 で節点 3 が閉じた）。

### ★★★★★★★★★★★★第 1166 で節点 3 が閉じた

`Found/GenEll/Lemma37Full.lean` の **`lemma_3_7_stableLine`**——
`Lemma 3.7` の第 3 の主張を **`HasLCyclicJ`（`Gal`-安定な直線）の側で**述べた形。

☆受けているのは `hdag`（`Lemma 3.5` の結論 `(†)` を安定直線の側で述べたもの）**だけ**であり、
第 1・第 2 の主張と第 3 の数値の核は `lemma_3_7`（第 1151）と**同じ材料で通った**。
★★すなわち **`Theorem 3.8` に繋げるのに残るのは `hdag` ただ 1 つ**であり、
それを埋めるのが節点 2d（埋め込み `L̄ ↪ L̄_v` の配管）である。

### ★★★★★★★★★★★★第 1154 で節点 1 が閉じた

`Found/GenEll/VeluDescent.lean` の **`veluQuotientFull_descends`**——
`W : WeierstrassCurve L` と `L̄` の座標の有限集合 `S` が**すべての `σ` で保たれる**なら、
`veluQuotientFull (W⁄L̄) S` は **`L` 上の曲線の底変換**である。

☆鍵は mathlib にそのものがあったことである——
`InfiniteGalois.mem_range_algebraMap_iff_fixed`（`L̄^{Gal} = L`）。
★標数 0 なら `IsGalois L L̄` はインスタンスで出る。

### ★★★★★★★★第 1153 で取れたもの

`Found/GenEll/VeluGaloisInvariant.lean`。★いずれも**無条件**である。

| 定理 | 内容 |
|---|---|
| `veluGx_semi` / `veluGy_semi` | 係数を固定する `Φ` に対し `Φ(veluGx(x,y)) = veluGx(Φx,Φy)` |
| `image_semiPair_eq` | `Φ` で保たれる有限集合は `Φ` の像で自分自身（単射だから） |
| `veluVFull_semi` / `veluWFull_semi` | 座標の集合が `Φ` で保たれるなら和は `Φ` で不変 |
| `fixesCoeffs_veluQuotientFull` | ★**商曲線の係数はすべて `Φ` で固定される** |

☆`/2` は `Φ 2 = 2` なので害が無い。★商の `a₁・a₂・a₃` は元のもののまま、
`a₄ = a₄ − 5v`・`a₆ = a₆ − b₂v − 7w` は `v`・`w` の不変性から出る。

### ★★★★節点 1 の道（測定）

`veluQuotientFull` は体 `F` と `S : Finset (F × F)` を受ける。★`F = L̄` で取れば
安定な直線の点の座標は `L̄` にあり、商 `veluCurve E v w` は `L̄` 上で定義される。
☆`H` が `Gal(L̄/L)`-安定なら `S` は `Gal` の作用で保たれるので、

    `v = Σ_{Q ∈ S} veluV2(Q)`,   `2w = Σ_{Q ∈ S} (veluU(Q) + 2 veluV2(Q)·x_Q)`

はいずれも `Gal` 不変、すなわち `L` の元である（`L̄^{Gal} = L`）。
★したがって `veluCurve E v w` は `L` 上の曲線の底変換である。

☆在庫: `Found/GenEll/VeluPointSet.lean` の `veluQuotientFull_image_eq`（座標の像で取った商）と
`Found/GaloisRep/VeluTateMuK.lean` の `veluQuotientFull_points_eq_field` が素材である。

### ★★★★★★★★★★★★節点 2 の道（第 1155 の測定）——**付値の議論が丸ごと要らない**

`Found/GaloisRep/TorsionIntegralGood.lean` は捻れ点の座標が
**`L` の付値環 `primeSubring p` に属する**ことを付値の言葉で示している
（`v(x) < 0` なら深さ `m` が取れて… という議論）。
☆安定直線の側では点の座標は `L̄` にあり、`L̄` に `p` の付値は一意には伸びない。

★★**だが `Lemma 3.5` が実際に要るのは Vélu の和 `v`・`w` が整であることだけ**であり、
第 1154 で**それらは `L` の元**だと分かった。☆したがって:

| 段 | 内容 | 状態 |
|---|---|---|
| 2a | `L̄` の捻れ点の `x` は `primeSubring p` 上**整** | ★第 1155（`isIntegral_x_of_addOrderOf_prime`） |
| 2b | `x` が整なら `y` も整（Weierstrass 方程式が `y` についてモニック） | ★第 1156（`isIntegral_y_of_isIntegral_x`） |
| 2c | `v` は多項式なので整、かつ `L` の元 ⟹ **属する** | ★第 1157（`isIntegral_veluVFull`・`mem_primeSubring_of_isIntegral_image`） |
| 2c' | `w` の `/2` | ★第 1158（`isIntegral_veluWFull_of_addOrderOf_prime`）——再添字して既存の対を使った |
| 2d | 悪い素点の側（`Δ_min` の関係）を安定直線で回す | ★**第 1179-1180 で `hcurveEq` の形まで取れた**（`veluQuotientFull_descends_algClosed`）。☆残るのは `ι(S)` が `μ_l` の座標集合になること（`Lemma 3.2, (i)`）を当てる段 |

★**付値の議論（深さ `m`・`ValAtLeast` の連鎖）がまるごと要らなくなる**のが本計測の利きである。

### ★★★★★★★★`w` の `/2` の道（第 1157 の測定）

☆`w` だけは割り算があるので多項式の議論では済まない。
★しかし**再添字すれば既存の道具がそのまま使える**:

`l`-巡回部分群 `H` は `L̄` の中で巡回であるから、生成元 `Q₀ ∈ E(L̄)` が取れて

    `S = ((range l).erase 0).image (fun k => pointCoords (k • Q₀))`

と書ける（`Q₀` は `L` 有理でなくてよい）。★反転は `k ↦ l − k` に対応するので、
`exists_veluW_of_inv`（第 960）と `exists_veluW_two`（第 1149）が**そのまま効く**。

☆つまり `TorsionIntegralGood.lean` の `isIntegral_veluQuotientFull_of_addOrderOf_prime`
を `L̄` と `A ≔ integralClosure R L̄` で組み直すだけである。
★**新しい数学は要らない**——Finset の対を作る道具を新規に書く必要は無い。

### ★★★★★★★★節点 2d（悪い素点の側）の実測（第 1159）

☆`Q` の**有理性**が悪い素点の連鎖で使われているのは、実は **1 か所だけ**である。

    `minDeltaExp_eq_mul_of_globalVelu'_K` の `hcurveEq`
    ——`E' ⊗ Lv` が Tate 曲線の `μ_l` による Vélu の商に一致すること

★局所の議論自体（`TateIsogenyK.lean`）は `ζ` を**完備体の中で**取っており、
`Q` が `L` 有理であることを**一度も使っていない**。
☆`hcurveEq` は `k • Q` の座標の `Lv` への像が `tatePhi(ζ^i)` と一致すること
（`Lemma 3.2, (i)`、`l ∤ v(q)` から出る）である。

★★**安定直線でも同じことが成り立つ**——`H` の点は `L̄` にあり、
`Lv` の代数閉包へ埋めれば同じ `μ_l` になる。
☆したがって 2d の中身は**埋め込み `L̄ ↪ L̄_v` の配管**であり、
数学的な新規は無い。★重みは配管の分だけである。

### ★★★★★★★★★★★★節点 2d の見積もりの**訂正**（第 1181）

★第 1159 で「`Q` の有理性が使われているのは `hcurveEq` の 1 か所だけ」と書き、
第 1179-1180 でその `hcurveEq` の形まで取った。☆**だが見積もりが甘かった**。

`TateIsogenyK.lean` の連鎖を読み直すと、`hcurveEq` の手前で

    `exists_primitiveRoot_of_torsion_point` → `ζ : p.adicCompletionIntegers L`

となっており、★★**`ζ` が `L_v` の中にあることを使っている**。
これは局所の点 `P` が `L_v` 有理（大域の `Q` が `L` 有理だから）であることに依る。

☆**安定直線の側では `ζ ∉ L_v` であり得る**ので、局所の議論を
`L_v(ζ_l)` へ上げ、そこで得た `Δ_min` の関係を `L_v` へ降ろす必要がある。

| 段 | 内容 | 素材 |
|---|---|---|
| 2d-1 | 局所の議論を `L_v(ζ_l)` で回す | 既存の `_K` の連鎖そのもの |
| 2d-2 | `Δ_min` の関係を `L_v` へ降ろす | ★★**第 1183 で閉じた**（`minDeltaExp_descend_of_baseChange`）。☆分岐指数が `l` と素であることは**要らなかった**——両辺とも `e` 倍になるので `e ≠ 0` で割るだけ |

★★★**素材はすでに `Found/GenEll/Thm38Kummer.lean` にある**——
それらはもともと `l ∤ v_K(q)` を `K(ζ_l)` へ上げるために書かれたものであり、
同じ分岐の帳簿がここでも使える。

☆重みの訂正: **4 → 8**。★新しい数学はないが、「1 か所の配管」ではない。

★★★★**再訂正（第 1183）**——**2d-2 は閉じた**。
`minDeltaExp_descend_of_baseChange`（`Found/GaloisRep/DegInfBaseChange.lean`）。
☆分岐指数が `l` と素であることは**要らなかった**——
上で `v_P(Δ_min(E')) = l·v_P(Δ_min(E))` なら、両辺とも `e` 倍なので
`e ≠ 0` で割れば下の関係が出る（`minDeltaExp_baseChange_of_semistableAt`、第 740）。
★重み 8 → **5**（残るのは 2d-1、局所の議論を `L_v(ζ_l)` で回す段）。

★★第 1185 で 2d-2 は**仮説なし**になった（`minDeltaExp_descend`）
——`e ≠ 0` は mathlib の `ramificationIdx_ne_zero_of_liesOver` から出る。

★第 1184 で 2d-1 の**入口**も取れた（`dvd_exponent_of_not_dvd_val`）
——`l ∤ v(q)` なら安定直線は `μ_l` である。
☆残るのはそれを受けて既存の `_K` の連鎖を `L_v(ζ_l)` で走らせる組み立てである。

#### ★★★★第 1186 の実測——**在庫に無い葉が 1 つある**

`_K` の連鎖（`minDeltaExp_eq_mul_at_bad_prime_K` 等）は
`hssE : SemistableAt p E`・`hssE' : SemistableAt p E'` を受ける。
★`L'` の上で走らせるには **`SemistableAt P (E.baseChange L')`** が要るが、
☆**それは在庫に無い**（`grep` で確認。底変換の補題は
`minDeltaExp_baseChange_of_semistableAt`・`minDeltaExp_eq_maxJ_baseChange` だけで、
どちらも下の半安定性を**仮説に取っている**）。

★道（測定済み）:

| 場合 | 道 |
|---|---|
| `minDeltaExp p E = 0` | `minDeltaExp_baseChange_le` と `minDeltaExp ≥ 0` で `= 0` |
| 極小モデル `C` で `v_p(c₄) = 0` | 同じ `C` を上げれば `v_P(c₄) = e·v_p(c₄) = 0`、`isMinimal_of_c4_vAdd_eq_zero`（第 320）で極小 |

☆分岐する場合の `valAdd` のスケーリングが要る
（`vAdd_algebraMap_eq_valAdd` は**不分岐を仮定している**、第 320）。
★`LocalHeightRamified.lean` の `ordAt_liesOver` がその役をする。

★★★★**第 1188 でその紐付けが取れた**——`valAdd_algebraMap_liesOver`
（`Found/GaloisRep/DegInfBaseChange.lean`）。
☆`valAdd` と `ordAt` は**同じ定義**（`-toAdd (unzero …)`）なので `rfl` で繋がり、
`ordAt_liesOver` がそのまま `valAdd P (algebraMap x) = e·valAdd p x` を与える。

★良い還元の場合は第 1187（`minDeltaExp_baseChange_eq_zero`）で取れているので、
☆残るのは極小モデルを上げて極小性を出す段である。

##### ★★★★第 1189 の実測——最後の葉は `primeSubring` 版の極小性判定

`isMinimal_of_c4_vAdd_eq_zero`（第 320、`Found/GaloisRep/LocalHeightDelta.lean:88`）は
**完備 DVR `R` と `tateDvrVal R K` の言葉**で書かれており、
`SemistableAt` が使う `IsMinimal (primeSubring p)` と `valAdd p` の言葉では**ない**。

★★**第 1190 でその版を取った**——`isMinimal_of_c4_valAdd_eq_zero`
（`Found/GaloisRep/LocalHeightDelta.lean`）。
☆`primeSubring p` は DVR で分数体が `L` なので第 320 がそのまま当てはまり、
橋は「2 つの加法付値が同時に `0` になる」だけである
（`valuation_isEquiv` ＋ mathlib の `Valuation.isEquiv_iff_val_eq_one`）。

★★★★**第 1192 で `SemistableAt` の底変換が閉じた**
——`semistableAt_baseChange`（`Found/GaloisRep/DegInfBaseChange.lean`）。
☆これで**節点 2d-1 の葉はすべて埋まった**。

##### ★★★★次の 1 本の名前（第 1193 の実測）

`minDeltaExp_eq_mul_of_globalVelu'_K` が `hQ`・`hE'` を使うのは
**`exists_point_j_tateModel`（`Found/GaloisRep/TateModelPoint.lean:58`）を呼ぶ 1 か所だけ**である。
その結論は「局所の点 `P` があって `j(E' ⊗ L_v)` が `P` による Vélu の商の `j` に等しい」である。

★★安定直線の側ではそれが**第 1180 ＋ 第 1184** から出る:
第 1180 で `E' ⊗ L̄_v = veluQuotientFull (E ⊗ L̄_v) (ι S)`、
第 1184 で `ι S` は `μ_l` の座標集合。☆`j` は底変換で保たれる。

##### ★★★★★★★★道の訂正（第 1194）——**`L(H)` へ上げて降ろす**

☆`exists_point_j_tateModel` の**結論は `P` が `L_v` 有理であることを要求する**。
安定直線の `μ_l` の点は `L_v(ζ_l)` 有理であって `L_v` 有理とは限らないので、
「`L_v(ζ_l)` で連鎖を走らせる」という局所の道は**そのままでは取れない**。

★★★**代わりに大域で `L'' ≔ L(H)`（`H` の点の定義体、次数は `l−1` を割る）へ上げる**:

| 段 | 中身 | 道具 |
|---|---|---|
| 1 | `L''` の上では `H` の生成元 `Q` が**有理** | 定義そのもの |
| 2 | `E' ⊗ L'' = veluQuotientFull (E ⊗ L'') (image (k•Q))` | ★★第 1197 `veluQuotientFull_descends_to_intermediate` |
| 3 | `E ⊗ L''`・`E' ⊗ L''` は半安定 | ★第 1192 `semistableAt_baseChange` |
| 4 | `L''` の上で既存の `_K` の連鎖が**そのまま走る** | 在庫 |
| 5 | 得た `Δ_min` の関係を `L` へ降ろす | ★第 1185 `minDeltaExp_descend` |

☆第 1152 で `L(H)` の道を「`ht^Falt`・`deg∞` が体を変えると動く」として退けたが、
★**その懸念は当たらない**——体を変えるのは**局所の `Δ_min` の計算だけ**であり、
大域の高さは `L` のままである（第 1185 で降りる）。

★★★★**第 1197 で 5 段すべての部品が揃った**——
段 1 第 1195、段 2 第 1197、段 3 第 1192、段 4 在庫、段 5 第 1185。
☆残るのはこれらを 1 本の `lemma_3_5_height_ineq_stableLine` に繋ぐ配管だけである。

##### ★★★★良い素点の側（`hfin`）も同じ道で降りる（第 1198 の実測）

☆`hfin` は**大域の高さの不等式**なので、体を変えると値が変わる
——「`L''` で走らせて降ろす」という道が使えないのではという懸念があった。

★**その懸念も当たらない**。在庫に 2 本ある:

| 成分 | 底変換での振る舞い | 在庫 |
|---|---|---|
| `∑ᶠ_p neronExp_p · log N(p)` | `[L'':L]` 倍 | ★`finsum_scaling`（`PrimeUnder.lean:232`、第 736） |
| `ht^Falt` | **変わらない**（半安定なら） | ★`htFaltOf_baseChange_of_semistable`（第 738） |

☆したがって `L''` で得た `hfin` を `[L'':L]` で割れば `L` の `hfin` になる。
★★**すなわち `L(H)` の道は `Lemma 3.5` 全体を覚える**——
悪い素点の側だけではない。

##### ★★★★★★★★★★★★第 1199——**降ろす段さえ要らなかった**

`Found/GaloisRep/Lemma35StableLine.lean` の
**`lemma_3_5_height_ineq_over_extension`**。

☆`Lemma 3.5` の不等式に現れる `ht^Falt` と `deg∞` は
**半安定な曲線では底変換で変わらない**（第 738・第 737）ので、
`L''` で得た不等式が**そのまま `L` の不等式**である。

★★★したがって `L''` の上で `Q` が有理でありさえすれば、
`L` の上の `Lemma 3.5` が出る——**節点 2 はこれで実質閉じた**。
☆残るのは `HasLCyclicJ`（`galRep` の安定直線）から
`L(H)` 上の有理な生成元を取り出す橋だけである。

##### ★★★★★★★★★★★★第 1201-1202——点と `T_l E` の第 1 層を結ぶ

| 定理 | 内容 | 第 |
|---|---|---|
| `tateProj_galTate` | ☆射影は `Gal` と可換（在庫、定義から `rfl`） | 205 |
| `addOrderOf_tateProj_one` | ★第 1 層の射影が `≠ 0` なら位数はちょうど `l` | 1201 |
| `exists_point_addOrderOf_eq_prime` | ★★★`L̄` 上に位数 `l` の点は**必ずある** | 1202 |
| `exists_tateProj_one_eq` | ★★★`l`-捽れ点は第 1 層に**持ち上がる** | 1202 |

★★★これで `Lemma 3.5` を有限拡大の上で回す（第 1199）ための
**生成元は必ず取れる**ことが言えた——`L̄` に取れば、
第 1195（`L(H)` は有限次）で有限拡大に落ちる。

☆残るのは `tateProj … 1` の**核**が `l • T_l E` であること
（`T_l E / l ≅ E[l]`）で、これが `HasLCyclicJ` の**直線**を
点の直線へ一意に対応させる段である。

☆（参考）部品の内訳——
第 1187（良い還元）・第 1188（分岐した `valAdd`）・第 1190（極小性）・
mathlib の `map_variableChange`。☆残るのは組み立てだけである。

### ☆別の道（採らない理由）

`L(H)`（`H` の点を有理化する体、次数は `l−1` を割る）へ上げる道もあるが、
★`ht^Falt`・`deg∞` が体を変えると動くので、`Lemma 3.7` の数値の帳簿を
書き直すことになる。☆**節点 1 の降下の方が安い**。

## ★★★★★★何が壊れていないか

☆第 1151 の `lemma_3_7`（`.src` は条なし）は**それ自体としては正しい**——
原文の第 1・第 2 の主張は完全に、第 3 の主張は (V) の読みで述べている。
★壊れているのは `Theorem 3.8` へ**繋ぐ**ところだけである。
★★`Lemma 3.5` も (V) の読みなので、そちらも節点 2 で言い直すことになる。
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta

/-! ## ★★★★★★★★片側の含意は立つ -/

/-- ★★★★★★★★**`L` 有理な生成元があれば `Gal`-安定な直線がある**（含意の向きの記録）。

☆逆は偽である——安定な直線の点は `L` 有理とは限らない。
★本ファイルはこの向きを**主張として証明しない**（`HasLCyclicJ` は `galRep` の
言葉で書かれており、`HasLCyclicVelu` の `Q` から直線を作るには
`Q` の `Tate` 加群への持ち上げが要る）。☆節点として記録するにとどめる。 -/
def lCyclicReadingFrame.src : Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(l-巡回の 2 つの読みが繋がっていない——Theorem 3.8 へ渡す枠)",
    sectionId := "genell-lemma-3-7" }

def lCyclicReadingFrame.needs : List ProofObligation :=
  [ .citation "[ABC3]" "lemma_3_7(第 1151、(V) の読みで証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.lemma_3_7") 1,
    .citation "[ABC3]" "imageContainsSL2J_of_alpha'(Theorem 3.8、残る仮説は α だけ)"
      (.inProject "ABC3" "ABC3.Found.GenEll.imageContainsSL2J_of_alpha'") 1,
    .implicitStep
      ("★★★★**2026-09-01（第 1152）の測定**——原文の `H_L` を " ++
       "(V)「`L` 有理点 `Q` で位数 `l`」と読んだのは `Lemma 3.5` が Vélu の商を作るためであり、" ++
       "`Theorem 3.8` が要るのは (J)「`Gal`-安定な直線」である。" ++
       "☆(V) ⟹ (J) は立つが逆は偽なので、`Lemma 3.7` の第 3 の主張の対偶が " ++
       "`Theorem 3.8` の `hno` に届かない。" ++
       "★直す道は節点 1（安定な直線による Vélu の商が `Gal` 不変性で `L` に降りること）→ " ++
       "節点 2（`Lemma 3.5` を安定直線の側で述べ直す）→ 節点 3（`Lemma 3.7` の言い直し）。") 26 ]

end ABC3.Skeleton.GenEll
