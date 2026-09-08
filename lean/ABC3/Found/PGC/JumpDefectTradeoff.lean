import ABC3.Found.PGC.EquivariantProjectionDescent

/-!
# [pGC] 跳び `t` と different `d` の**トレードオフは等式**である —— ただし `p ≥ 5` では届かない

配られた持ち場は

> ★★★`d_{M/E}` と跳び `t` のトレードオフ —— `AxLemmaGraded` に残った唯一の数学。
> これが出れば `towerBudget_le` に代入するだけで `AxLemmaGraded` → `AxSenTate` が閉じる。

であった。★★**結論を先に書く: トレードオフは出た(しかも等式である)。しかし
「代入するだけで閉じる」は `p ≥ 5` で偽である。** 以下、測定を先に書く。

## ★★★測定 0 —— トレードオフの正体は**等式**(見立ては半分外れ)

本体の見立ては「different の指数 `d = (p−1)(t+1)` と `‖P‖ = p^{(d−e+1)/e}` を
突き合わせるだけ」であった。★`‖P‖` の式は **`d` の符号が逆**である。正しくは
Serre `Corps Locaux` V §3 Lemme 4(`Tr_{M/E}(𝔭_M^n) = 𝔭_E^{⌊(n+d)/e⌋}`)から

  ★`‖P‖ ≤ p^{γ/e_M}`,  `γ := v_M(p) + (p−1) − d_{M/E}`   (`P = (1/p)Tr_{M/E}`)

であり、`d_{M/E} = (p−1)(t+1)`(Serre IV §1 Prop 4、`M/E` は `p` 次で跳びは 1 個)を入れると
**`(p−1)` がちょうど消えて**

  ★★★`γ = v_M(p) − (p−1)·t`   (**等式**。不等式ではない)

となる。★これがトレードオフの正体である。`RamificationJumpBound` の sharp な
`(p−1)i ≤ e_L` は「`γ ≥ 0`」と**同じ 1 つの不等式**であって、別の情報ではない。
直ちに

  ★★`t + γ = v_M(p) − (p−2)·t ≤ v_M(p)`   (等号は `p = 2` または `t = 0` のときだけ)

が出る(`JumpDefect.jump_add_defect_le`)。★`γ = p−1` が円分層 4 本で一致していた
(`EquivariantProjectionDescent` の測定 1)のも、この等式に `t = p^{n−1} − 1`,
`v_M(p) = p^{n−1}(p−1)` を入れれば `γ = p^{n−1}(p−1) − (p−1)(p^{n−1}−1) = p−1` と出る。

## ★★★測定 1 —— しかし `t + γ ≤ v_M(p)` は**予算に届かない**

対の予算は `v_M` の目盛りで `B_k = e_K·p(p^k−1)/(p−1)²` である
(`PairLedger.rpow_div_le_prod_Icc_axDecay_iff` に `E = e_M = p^k e_K` を入れた形)。
`k = 2`, `p = 3`, `e_K = 2` で `B_2 = 12` に対し `v_M(p) = e_M = 18` なので、
★**`t + γ ≤ v_M(p)` は予算より弱い。トレードオフ単独では閉じない。**

## ★★★測定 2 —— 2 層の閉じた形(★既存ファイルの数値をすべて再現した)

`k = 2` の塔 `K ⊂ F_1 ⊂ M`(各層は全分岐巡回 `p` 次)を、`Γ = Gal(M/K)` の
下付き跳び `S := t_1 < T := t_2` で書く(★`t_j` は同時に「層 `F_j/F_{j−1}` の跳びを
`v_{F_j}` で測った値」でもある)。`v_M` の目盛りで:

| 量 | 式 |
|---|---|
| 上の層 `M/F_1` の降下 | `T` |
| 下の層 `F_1/K` の降下 | `p·S` |
| 射影 `(1/p)Tr_{M/F_1}` の欠損 | `γ = p²e_K − (p−1)T` |
| `σ−1` の収縮(★全 `Γ` で有効) | `S` |
| ★収縮つき合成 | `C = T + (p−1)S` |
| ★射影つき合成 | `P = max(T, pS) + γ` |
| 対の予算 | `B = e_K·p(p+1)/(p−1)` |

★★安定域(下記 `hstab`)では `P = 2p·e_K − (p−2)·S` と**閉じた形**になる
(`JumpDefect.projCost_eq_of_stable`)。既存ファイルの実測値はすべてこれで出る:

| 塔 | `p` | `e_K` | `S` | `T` | `C` | `P` | `B` | 既存ファイルの記述 |
|---|---|---|---|---|---|---|---|---|
| `ℚ₃(ζ₂₇)/ℚ₃(ζ₉)/ℚ₃(ζ₃)` | 3 | 2 | 2 | 8 | ★`12` | ★`10` | `12` | 収縮 `12` / 射影 `10` / 予算 `12` ✓ |
| `ℚ₃(ζ₈₁)/ℚ₃(ζ₂₇)/ℚ₃(ζ₉)` | 3 | 6 | 8 | 26 | ★`42` | ★`28` | `36` | 収縮 `42` / 射影 `28` / 予算 `36` ✓ |

★`p = 3` では **slack がちょうど `S`(第 1 跳び)**になる(`JumpDefect.slack_three_eq`):
`12 − 10 = 2 = S`、`36 − 28 = 8 = S`。★偶然ではない。

## ★★★★測定 3 —— **`p ≥ 5` で 2 つの枝の間に窓が空く**(★本ファイルの主結果)

跳びの列に効く制約は 2 つだけである(★どちらも古典):

* **層の上界** `(p−1)S ≤ p·e_K`、`(p−1)T ≤ p²·e_K`(Serre III §6 Prop 13、`RamificationJumpBound`)
* ★**二分律**: `(p−1)S ≤ e_K`(臨界 `e/(p−1)` 以下) または `T = S + p·e_K`(安定域)。
  後者は `n > e/(p−1)` で `(U^{(n)})^p = U^{(n+e)}` から出る古典
  (上付き跳びは `u_2 = u_1 + e_K`、Herbrand で `T = S + p·e_K`)。

この 2 つだけから:

| 場合 | 使う枝 | 判定 |
|---|---|---|
| `(p−1)S ≤ e_K` | 収縮 `C` | ★**すべての `p` で入る** (`(p−1)C ≤ (p²+p−1)e < p(p+1)e`) |
| `T = S + p e_K` かつ `(p−1)S ≤ 2e_K` | 収縮 `C` | 入る |
| `T = S + p e_K` かつ `(p−1)(p−2)S ≥ p e_K(p−3)` | 射影 `P` | 入る |

★★2 つの十分条件が**全体を覆うのは `2(p−2) ≥ p(p−3)`、すなわち `(p−1)(p−4) ≤ 0`、
すなわち `p ≤ 3` のとき、かつそのときに限る**(`JumpDefect.covered_iff` / `covered_of_le_three`)。

⇒ ★★★**`p ∈ {2,3}`: この道は閉じる**(`JumpDefect.min_cost_fits_of_le_three`)。
⇒ ★★★**`p ≥ 5`: 窓 `2e/(p−1) < S < p e(p−3)/((p−1)(p−2))` が空で、どちらの枝も届かない。**

★最小の整数点は `p = 5`, `e_K = 3`, `S = 2`, `T = 17`(⇒ `C = 25`, `P = 24`, `B = 22.5`)。
★★★一般に **`(e_K, S, T) = (p−2, 2, 2 + p(p−2))` が `p ≥ 5` のすべてで反例になる**
(`JumpDefect.witness_gap_exceeds`。収縮 `C = p²`、射影 `P = 2(p−2)(p−1)`、
予算の台帳 `p(p+1)(p−2)` に対して `p³−p² > p³−p²−2p` と `p³−7p²+12p−4 > 0`)。
★★この `(S,T)` は許容である(層の上界 `8 ≤ 15`, `68 ≤ 75`、安定域 `17 = 2 + 15` を満たす)。
`JumpDefect.witness_five_exceeds_budget` / `JumpDefect.witness_five_exceeds_pair_budget`。

★★★**しかもこの `24` は古典の Ax の定数 `p^{p/(p−1)²}`(= `v_M` で `375/16 = 23.4375`)
すら超える**(`JumpDefect.witness_five_exceeds_ax_constant`)。Ax–Sen–Tate は古典の定理なので、
⇒ ★★★**「偽なのは定理ではなく、この上界の取り方である」**。`p ≥ 5` では
「sharp な層の限界 ＋ 同変射影 ＋ 収縮」の 3 つを最善に組み合わせても足りない。

★総当たり(厳密有理数、`.../scratchpad/trade/trade.py`, `trade2.py`, `closed2.py`):
`p ≤ 3`, `e ≤ 24`, `k ≤ 4` の許容列 43,188 本で ★**破れ 0 件**、
`p = 5,7,11,13` では `k = 2` に限っても破れが出る(`p=5, e≤30` で 118 件)。

## 本ファイルが出したもの(すべて `sorry` 0、`Found/` に `sorry` を残さない)

| 宣言 | 内容 |
|---|---|
| ★★★`JumpDefect.defect_eq` | 抽象核: 欠損は `V − (p−1)t`(★**等式**。トレードオフの正体) |
| ★★`JumpDefect.jump_add_defect_le` / `_eq_iff` | `t + γ ≤ V`、等号は `p = 2 ∨ t = 0` |
| ★★`JumpDefect.contrCost_fits` | 収縮の枝は `(p−1)S ≤ e` で入る(★`p` に依らない) |
| ★★`JumpDefect.projCost_eq_of_stable` | 射影の枝の閉じた形 `2pe − (p−2)S` |
| ★★★`JumpDefect.min_cost_fits_of_le_three` | ★`p ≤ 3` なら 2 枝の `min` は予算に入る |
| ★★★`JumpDefect.gap_pos_of_five_le` / `covered_iff` | ★`p ≥ 5` で窓が空く(`(p−1)(p−4) > 0`) |
| ★★★`JumpDefect.witness_five_*` | ★`p=5, e=3, S=2, T=17` は許容だが予算も Ax の定数も超える |
| ★★★★`JumpDefect.witness_gap_exceeds` | ★**`p ≥ 5` のすべて**で反例の族 `(p−2, 2, 2+p(p−2))` |
| ★★`JumpDefect.pairBudget_iff` | 台帳: `(p−1)J ≤ p(p+1)e` ⇔ `p^{J/(p²e)} ≤ ∏_{i∈[1,2]} axDecay p i` |
| `JumpDefect.Zeta27.*` | 既存の実測値を閉じた形が再現する |

★抽象核 §1・§2 は **`ℤ` の算術だけ**で、分岐・付値・Galois の語彙が 1 語も出ない。

## ★★★埋まっていないもの(★正確に)

1. ★**`p ≥ 5` の `AxLemmaGraded` は本ファイルの道では閉じない**。上の窓の中の塔で、
   sharp な層の限界 + 同変射影 + 収縮の**どの組み合わせでも**対の予算を超える。
   ★閉じるには別の機構(層の限界そのものが sharp でない、という改良)が要る。
   ★★測定が示唆する正しい形(★**予想。本ファイルは証明していない**): 循環な塔では
   実測の最悪値が**上の層の跳び `T` そのもの**になっている。
   ★確かめられたのは 2 本だけである(既存ファイルが `L` を明示しているもの):
   `ℚ₃(ζ₂₇)/ℚ₃(ζ₃)` の `L = 8 = T`、`ℚ₃(ζ₈₁)/ℚ₃(ζ₉)` の `L = 26 = T`。
   ★`ℚ₂` 系は `DeepDescentPairDirect` の表が slack しか書いておらず、対の `d` が深さと
   一致しているかが読み取れないので**数えていない**。
   ★★逆に、**循環でない塔ではこの予想は危うい**: `ℚ₂(ζ₁₆)/ℚ₂`(`Gal ≅ ℤ/2×ℤ/4`)は
   最大の跳びが `7` だが、同表の slack から逆算すると `8` になる(★この逆算自体が
   上の理由で確かめられていない)。★予想を使う前にここを測ること。
   ★なお `T ≤ B_k` は `k ≥ 2` なら常に真である(`(p/(p−1))(1−p^{−k}) > 1`)ので、
   予想が真なら `AxLemmaGraded` は循環な塔では従う。
2. `p ≤ 3` でも、本ファイルが閉じたのは `k = 2`(対)だけである。`k ≥ 3` は総当たりで
   破れ 0 件を測ったが、閉じた形の証明は書いていない。
3. 二分律(安定域)は仮説 `hstab` として置いた。木にも mathlib にも無い
   (`grep -n "U\^{(n)}\|stable range" lean/ABC3/Found/PGC/*.lean` で 0 件)。

## 逸脱の記録(CLAUDE.md「逸脱」)

1. 抽象核は `ℤ` の不等式で、`p` の素数性を使わない(`2 ≤ p` のみ)。分岐の側から来る
   仮説(層の上界・二分律)は**仮定として置く**。★これは `EquivariantProjectionDescent` の
   逸脱 2(測定は Lean の外)と同じ扱いである。
2. `γ = v_M(p) − (p−1)t` の導出に Serre V §3 Lemme 4 と IV §1 Prop 4 を使うが、
   ★どちらも木にも mathlib にも無い(`differentIdeal` はあるが跳びとの関係は無い)。
   本ファイルはこの 2 つを**式の形で仮定に取り込んだ**(`defect_eq` の仮説 `hd`)。
3. 収縮の枝 `C = T + (p−1)S` は「`σ−1` の収縮は全 `Γ` で `S`」を使う。深い `σ` では
   もっと縮む(`T`)が、★下の層の降下に要るのは**浅い `σ`** なので `S` が正しい。
-/

namespace ABC3.Found.PGC

/-! ## §1 抽象核 —— トレードオフは**等式**である(`ℤ` の算術のみ) -/

namespace JumpDefect

/-- ★★★**抽象核 1(トレードオフの正体)**。分岐・付値・Galois の語彙が 1 語も出ない。

`d = (p−1)(t+1)` のとき、欠損 `V + (p−1) − d` は ★**ちょうど `V − (p−1)t`** になる。
`(p−1)` が消えるのがこの補題の全部である。

具体層での読み: `V = v_M(p)`、`d = d_{M/E}`(different の指数)、`t` = `M/E` の跳び。
左辺は `‖(1/p)Tr_{M/E}‖ = p^{(V + (p−1) − d)/e_M}` の肩(Serre V §3 Lemme 4)、
`d = (p−1)(t+1)` は Serre IV §1 Prop 4。 -/
theorem defect_eq (p t V d : ℤ) (hd : d = (p - 1) * (t + 1)) :
    V + (p - 1) - d = V - (p - 1) * t := by
  subst hd; ring

/-- ★★`t + γ = V − (p−2)t` —— 跳びと欠損の**和**の閉じた形。 -/
theorem jump_add_defect (p t V : ℤ) : t + (V - (p - 1) * t) = V - (p - 2) * t := by ring

/-- ★★★**トレードオフ**(配られた持ち場そのもの): `t + γ ≤ V`。
★`t` が小さいと `γ` が大きく、`t` が大きいと `γ` が小さい —— その和は `V` を超えない。 -/
theorem jump_add_defect_le {p t V : ℤ} (hp : 2 ≤ p) (ht : 0 ≤ t) :
    t + (V - (p - 1) * t) ≤ V := by nlinarith

/-- ★等号は `p = 2` か `t = 0` のときだけ(★`p ≥ 3` では真に得をする)。 -/
theorem jump_add_defect_eq_iff {p t V : ℤ} (hp : 2 ≤ p) (ht : 0 ≤ t) :
    t + (V - (p - 1) * t) = V ↔ (p = 2 ∨ t = 0) := by
  constructor
  · intro h
    have h2 : (p - 2) * t = 0 := by linarith [jump_add_defect p t V]
    rcases mul_eq_zero.mp h2 with h3 | h3
    · exact Or.inl (by linarith)
    · exact Or.inr h3
  · rintro (rfl | rfl) <;> ring

/-- ★`γ ≥ 0` は sharp な跳びの上界 `(p−1)t ≤ V`(`RamificationJumpBound`)と**同値**であって、
別の情報ではない。 -/
theorem defect_nonneg_iff {p t V : ℤ} : 0 ≤ V - (p - 1) * t ↔ (p - 1) * t ≤ V := by
  constructor <;> intro h <;> linarith

/-! ## §2 抽象核 2 —— 2 つの枝と対の予算(★`ℤ` の算術のみ、分岐の語彙なし) -/

/-- ★**収縮の枝**。上の層を `T` 払って降り、`σ−1` の収縮 `S` の分だけ `ε` の膨らみが減り、
下の層を `p·S` 払う: `max(T, pS + (T − S)) = T + (p−1)S`(`S ≤ T` のとき)。 -/
def contrCost (p S T : ℤ) : ℤ := T + (p - 1) * S

/-- ★**射影の枝**。同変射影の欠損 `γ = p²e − (p−1)T`(§1 の `defect_eq` に `V = p²e` を入れた形)を
`EquivProj.exists_of_proj_two_descent` の `max (A·max 1 Q) (B·Q)` に入れた指数。 -/
def projCost (p e S T : ℤ) : ℤ := max T (p * S) + (p ^ 2 * e - (p - 1) * T)

/-- ★**対の予算を `(p−1)` 倍して分母を払った形**。
`PairLedger.rpow_div_le_prod_Icc_axDecay_iff` に `E = e_M = p²e`, `k' = 0`, `d = 2` を入れると
`(p−1)²p²·J ≤ (p²−1)p·p²e`、すなわち `(p−1)·J ≤ p(p+1)e`。 -/
def pairBudget (p e : ℤ) : ℤ := p * (p + 1) * e

/-- ★★**収縮の枝は臨界以下(`(p−1)S ≤ e`)なら必ず入る** —— ★`p` に依らない。
`(p−1)C = (p−1)T + (p−1)²S ≤ p²e + (p−1)e = (p²+p−1)e < p(p+1)e`。 -/
theorem contrCost_fits {p e S T : ℤ} (hp : 2 ≤ p) (he : 0 ≤ e)
    (hS : (p - 1) * S ≤ e) (hT : (p - 1) * T ≤ p ^ 2 * e) :
    (p - 1) * contrCost p S T ≤ pairBudget p e := by
  have h1 : (0:ℤ) ≤ p - 1 := by linarith
  have h2 : (p - 1) * ((p - 1) * S) ≤ (p - 1) * e := mul_le_mul_of_nonneg_left hS h1
  unfold contrCost pairBudget
  nlinarith

/-- ★★**射影の枝の閉じた形**(安定域 `T = S + pe`): `P = 2pe − (p−2)S`。
★これで `ℚ₃(ζ₂₇)` の `10`、`ℚ₃(ζ₈₁)` の `28` がそのまま出る。 -/
theorem projCost_eq_of_stable {p e S T : ℤ} (hstab : T = S + p * e)
    (hSle : (p - 1) * S ≤ p * e) :
    projCost p e S T = 2 * p * e - (p - 2) * S := by
  subst hstab
  unfold projCost
  rw [max_eq_left (by nlinarith : p * S ≤ S + p * e)]
  ring

/-- ★★**安定域でも収縮の枝が入る条件**は `(p−1)S ≤ 2e`(`C = pS + pe`)。 -/
theorem contrCost_fits_of_stable {p e S T : ℤ} (hp : 2 ≤ p)
    (hstab : T = S + p * e) (hS2 : (p - 1) * S ≤ 2 * e) :
    (p - 1) * contrCost p S T ≤ pairBudget p e := by
  subst hstab
  have h1 : (0:ℤ) ≤ p := by linarith
  have h2 : p * ((p - 1) * S) ≤ p * (2 * e) := mul_le_mul_of_nonneg_left hS2 h1
  unfold contrCost pairBudget
  nlinarith

/-- ★★**安定域で射影の枝が入る条件**は `p e (p−3) ≤ (p−1)(p−2)S`。 -/
theorem projCost_fits_iff_of_stable {p e S T : ℤ}
    (hstab : T = S + p * e) (hSle : (p - 1) * S ≤ p * e) :
    (p - 1) * projCost p e S T ≤ pairBudget p e ↔ p * e * (p - 3) ≤ (p - 1) * ((p - 2) * S) := by
  rw [projCost_eq_of_stable hstab hSle]
  unfold pairBudget
  constructor <;> intro h <;> nlinarith

/-- ★★★**`p ≤ 3` なら安定域でも射影の枝が入る**(`pe(p−3) ≤ 0 ≤ (p−1)(p−2)S`)。 -/
theorem projCost_fits_of_le_three {p e S T : ℤ} (hp : 2 ≤ p) (hp3 : p ≤ 3)
    (hS : 0 ≤ S) (hstab : T = S + p * e) (hSle : (p - 1) * S ≤ p * e) :
    (p - 1) * projCost p e S T ≤ pairBudget p e := by
  rw [projCost_fits_iff_of_stable hstab hSle]
  have h1 : (0:ℤ) ≤ (p - 1) * ((p - 2) * S) :=
    mul_nonneg (by linarith) (mul_nonneg (by linarith) hS)
  nlinarith

/-- ★★★★**本ファイルの主結果(肯定側)**: `p ≤ 3` なら、層の上界と二分律だけから
2 つの枝の `min` は**対の予算に入る**。

仮説はすべて古典:
* `hSle` / `hTle` —— 層ごとの sharp な跳びの上界(`RamificationJumpBound` の `(p−1)i ≤ e_L`)
* `hdich` —— ★二分律「臨界 `e/(p−1)` 以下」または「安定域 `T = S + pe`」

★`p ≥ 5` ではこの結論は**偽**である(`witness_five_exceeds_budget`)。 -/
theorem min_cost_fits_of_le_three {p e S T : ℤ} (hp : 2 ≤ p) (hp3 : p ≤ 3) (he : 0 ≤ e)
    (hS : 0 ≤ S) (hSle : (p - 1) * S ≤ p * e) (hTle : (p - 1) * T ≤ p ^ 2 * e)
    (hdich : (p - 1) * S ≤ e ∨ T = S + p * e) :
    (p - 1) * min (contrCost p S T) (projCost p e S T) ≤ pairBudget p e := by
  have h1 : (0:ℤ) ≤ p - 1 := by linarith
  rcases hdich with h | h
  · exact le_trans (mul_le_mul_of_nonneg_left (min_le_left _ _) h1)
      (contrCost_fits hp he h hTle)
  · exact le_trans (mul_le_mul_of_nonneg_left (min_le_right _ _) h1)
      (projCost_fits_of_le_three hp hp3 hS h hSle)

/-! ### ★§2.1 2 つの十分条件が全体を覆うのは `p ≤ 3` のときだけ -/

/-- ★★**覆う条件は `(p−1)(p−4) ≤ 0`**。左辺は「射影が要求する `S` の下限」、
右辺は「収縮が許す `S` の上限」を `e/((p−1)(p−2))` 倍したもの。 -/
theorem covered_iff {p : ℤ} : p * (p - 3) ≤ 2 * (p - 2) ↔ (p - 1) * (p - 4) ≤ 0 := by
  constructor <;> intro h <;> nlinarith

/-- ★`p ≤ 3` なら覆う。 -/
theorem covered_of_le_three {p : ℤ} (hp : 2 ≤ p) (hp3 : p ≤ 3) : p * (p - 3) ≤ 2 * (p - 2) := by
  nlinarith

/-- ★★★**`p ≥ 5` では窓が空く**(`0 < (p−1)(p−4)`)。
すなわち `2e/(p−1) < S < p e (p−3)/((p−1)(p−2))` を満たす `S` が存在し、
そこでは収縮も射影も予算に届かない。 -/
theorem gap_pos_of_five_le {p : ℤ} (hp : 5 ≤ p) : 2 * (p - 2) < p * (p - 3) := by nlinarith

/-! ## §3 ★★★反例(`p = 5`) —— 許容な塔で**どちらの枝も予算を超える** -/

/-- ★`p = 5`, `e_K = 3`, `S = 2`, `T = 17` は**許容**である:
層の上界 `(p−1)S = 8 ≤ 15 = pe`、`(p−1)T = 68 ≤ 75 = p²e`、
二分律は安定域の側 `T = S + pe = 2 + 15 = 17`(臨界の側 `(p−1)S ≤ e` は成り立たない)。 -/
theorem witness_five_admissible :
    (5 - 1) * (2:ℤ) ≤ 5 * 3 ∧ (5 - 1) * (17:ℤ) ≤ 5 ^ 2 * 3 ∧
      (17:ℤ) = 2 + 5 * 3 ∧ ¬ ((5 - 1) * (2:ℤ) ≤ 3) := by norm_num

/-- ★2 つの枝の値: 収縮 `25`、射影 `24`。 -/
theorem witness_five_costs :
    contrCost 5 2 17 = 25 ∧ projCost 5 3 2 17 = 24 ∧ pairBudget 5 3 = 90 := by
  refine ⟨by norm_num [contrCost], ?_, by norm_num [pairBudget]⟩
  norm_num [projCost]

/-- ★★★**主結果(否定側)**: `p = 5` の許容な塔で `(p−1)·min(C,P) = 96 > 90 = ` 予算。
⇒ ★「トレードオフを `towerBudget_le` に代入するだけで `AxLemmaGraded` が閉じる」は**偽**。 -/
theorem witness_five_exceeds_budget :
    pairBudget 5 3 < (5 - 1) * min (contrCost 5 2 17) (projCost 5 3 2 17) := by
  norm_num [pairBudget, contrCost, projCost]

/-- ★★★**しかも古典の Ax の定数 `p^{p/(p−1)²}` すら超える**:
`v_M` の目盛りで Ax の上界は `p·v_M(p)/(p−1)² = 5·75/16 = 375/16 = 23.4375` だが、
この道の上界は `24` である(`(p−1)²·24 = 384 > 375 = p·p²e`)。
★Ax–Sen–Tate は古典の定理なので、⇒ ★**偽なのは定理ではなくこの上界の取り方**である。 -/
theorem witness_five_exceeds_ax_constant :
    5 * (5 ^ 2 * 3) < (5 - 1) ^ 2 * min (contrCost 5 2 17) (projCost 5 3 2 17) := by
  norm_num [contrCost, projCost]

/-! ### ★§3.1 反例は `p = 5` 限りではない —— **`p ≥ 5` のすべてで族が取れる** -/

/-- ★`e_K = p−2`, `S = 2`, `T = 2 + p(p−2)` は**許容**である(`p ≥ 5`):
層の上界 2 つと、二分律の「安定域」側を満たし、「臨界以下」側は満たさない。 -/
theorem witness_gap_admissible {p : ℤ} (hp : 5 ≤ p) :
    (p - 1) * 2 ≤ p * (p - 2) ∧ (p - 1) * (2 + p * (p - 2)) ≤ p ^ 2 * (p - 2) ∧
      ¬ ((p - 1) * 2 ≤ p - 2) := by
  refine ⟨by nlinarith, by nlinarith, ?_⟩
  intro h; linarith

/-- ★収縮の枝は `p²` になる。 -/
theorem witness_gap_contrCost (p : ℤ) : contrCost p 2 (2 + p * (p - 2)) = p ^ 2 := by
  unfold contrCost; ring

/-- ★射影の枝は `2(p−2)(p−1)` になる。 -/
theorem witness_gap_projCost {p : ℤ} (hp : 5 ≤ p) :
    projCost p (p - 2) 2 (2 + p * (p - 2)) = 2 * (p - 2) * (p - 1) := by
  rw [projCost_eq_of_stable (by ring) (by nlinarith)]
  ring

/-- ★★★★**主結果(否定側・一般形)**: `p ≥ 5` なら、許容な `(e_K, S, T) = (p−2, 2, 2+p(p−2))` で
★**収縮も射影も対の予算を超える**(`p³−p² > p³−p²−2p` と `p³−7p²+12p−4 > 0`)。

⇒ ★★**`p ≥ 5` では「sharp な層の限界 + 同変射影 + 収縮」だけでは `AxLemmaGraded` は閉じない。**
`p = 5` で `e_K = 3`, `S = 2`, `T = 17`, `C = 25`, `P = 24`, 予算 `22.5`。 -/
theorem witness_gap_exceeds {p : ℤ} (hp : 5 ≤ p) :
    pairBudget p (p - 2)
      < (p - 1) * min (contrCost p 2 (2 + p * (p - 2))) (projCost p (p - 2) 2 (2 + p * (p - 2))) := by
  have h5 : (0:ℤ) ≤ p - 5 := by linarith
  rw [witness_gap_contrCost, witness_gap_projCost hp]
  have hC : pairBudget p (p - 2) < (p - 1) * p ^ 2 := by
    unfold pairBudget; nlinarith
  have hP : pairBudget p (p - 2) < (p - 1) * (2 * (p - 2) * (p - 1)) := by
    unfold pairBudget
    nlinarith [mul_nonneg (mul_nonneg h5 h5) h5, mul_nonneg h5 h5]
  rcases min_cases (p ^ 2) (2 * (p - 2) * (p - 1)) with ⟨h, _⟩ | ⟨h, _⟩ <;> rw [h] <;> assumption

/-! ## §4 ★台帳 —— `ℤ` の不等式と `axDecay` の積を繋ぐ -/

/-- ★★**対の予算の台帳**(`k' = 0`, `d = 2`, `E = e_M = p²e`)。
`PairLedger.rpow_div_le_prod_Icc_axDecay_iff` から `(p−1)p²` を約せば
★`(p−1)·J ≤ p(p+1)e` に簡約される。これが `pairBudget` の定義の根拠である。 -/
theorem pairBudget_iff {p : ℕ} [Fact p.Prime] {J e : ℕ} (he : 0 < e) :
    (p : ℝ) ^ ((J : ℝ) / ((p ^ 2 * e : ℕ) : ℝ)) ≤ ∏ i ∈ Finset.Icc 1 2, axDecay p i
      ↔ ((p : ℤ) - 1) * (J : ℤ) ≤ pairBudget (p : ℤ) (e : ℤ) := by
  have hp2 : 2 ≤ p := (Fact.out : p.Prime).two_le
  have hE : 0 < p ^ 2 * e := by positivity
  have hIcc : (Finset.Icc 1 2 : Finset ℕ) = Finset.Icc (0 + 1) (0 + 2) := by norm_num
  rw [hIcc, PairLedger.rpow_div_le_prod_Icc_axDecay_iff (p := p) (k' := 0) (d := 2)
    (J := J) (E := p ^ 2 * e) hE]
  have hp2' : (2:ℤ) ≤ (p:ℤ) := by exact_mod_cast hp2
  have h1 : 1 ≤ p := by omega
  have h2 : 1 ≤ p ^ 2 := Nat.one_le_pow _ _ (by omega)
  unfold pairBudget
  zify [h1, h2]
  have hpos : (0:ℤ) < ((p:ℤ) - 1) * (p:ℤ) ^ 2 :=
    mul_pos (by linarith) (by positivity)
  have e1 : ((p:ℤ) - 1) ^ 2 * (p:ℤ) ^ (0 + 2) * (J:ℤ)
      = (((p:ℤ) - 1) * (J:ℤ)) * (((p:ℤ) - 1) * (p:ℤ) ^ 2) := by ring
  have e2 : ((p:ℤ) ^ 2 - 1) * (p:ℤ) * ((p:ℤ) ^ 2 * (e:ℤ))
      = ((p:ℤ) * ((p:ℤ) + 1) * (e:ℤ)) * (((p:ℤ) - 1) * (p:ℤ) ^ 2) := by ring
  rw [e1, e2, mul_le_mul_iff_left₀ hpos]

/-! ## §5 数値層 —— 既存ファイルの実測値を閉じた形が再現する -/

namespace Zeta27

/-- ★`ℚ₃(ζ₂₇) ⊃ ℚ₃(ζ₉) ⊃ ℚ₃(ζ₃)`(`p = 3`, `e_K = 2`, `S = 2`, `T = 8`)の 3 つの値。
★`EquivariantProjectionDescent` の測定 2(収縮 `12` / 射影 `10` / 予算 `12`)と一致する。 -/
theorem costs : contrCost 3 2 8 = 12 ∧ projCost 3 2 2 8 = 10 ∧ pairBudget 3 2 = 2 * 12 := by
  refine ⟨by norm_num [contrCost], ?_, by norm_num [pairBudget]⟩
  norm_num [projCost]

/-- ★射影の枝は予算に入る(`10 ≤ 12`)。 -/
theorem proj_fits : (3 - 1) * projCost 3 2 2 8 ≤ pairBudget 3 2 := by
  norm_num [projCost, pairBudget]

/-- ★収縮の枝は**ちょうど**予算(`12 = 12`)—— `p = 3` の slack `S = 2` の分だけ射影が得をしている。 -/
theorem contr_eq_budget : (3 - 1) * contrCost 3 2 8 = pairBudget 3 2 := by
  norm_num [contrCost, pairBudget]

/-- ★`ℚ₃(ζ₈₁) ⊃ ℚ₃(ζ₂₇) ⊃ ℚ₃(ζ₉)`(`e_K = 6`, `S = 8`, `T = 26`): 収縮 `42` は
予算 `36` を**超える**が、射影 `28` は入る。★`EquivariantProjectionDescent` の測定 2 と一致。 -/
theorem zeta81_costs :
    contrCost 3 8 26 = 42 ∧ projCost 3 6 8 26 = 28 ∧ pairBudget 3 6 = 2 * 36 := by
  refine ⟨by norm_num [contrCost], ?_, by norm_num [pairBudget]⟩
  norm_num [projCost]

theorem zeta81_contr_exceeds : pairBudget 3 6 < (3 - 1) * contrCost 3 8 26 := by
  norm_num [contrCost, pairBudget]

theorem zeta81_proj_fits : (3 - 1) * projCost 3 6 8 26 ≤ pairBudget 3 6 := by
  norm_num [projCost, pairBudget]

end Zeta27

/-- ★★★**`p = 3` の slack はちょうど第 1 跳び `S`**(台帳の目盛りで `(p−1)S = 2S`)。
`ℚ₃(ζ₂₇)`: `24 − 2·10 = 4 = 2·2`、`ℚ₃(ζ₈₁)`: `72 − 2·28 = 16 = 2·8`。
★★これは `p = 3` の射影の枝が予算に**必ず**入ることの定量版である。 -/
theorem slack_three_eq {e S T : ℤ} (hstab : T = S + 3 * e) (hSle : 2 * S ≤ 3 * e) :
    pairBudget 3 e - (3 - 1) * projCost 3 e S T = 2 * S := by
  rw [projCost_eq_of_stable (p := 3) hstab (by linarith)]
  unfold pairBudget
  ring

/-- ★★★**反例の rpow 版**: `p = 5`, `e_M = 75`, `J = 24` は対の予算に**入らない**
(`(p−1)²p²J = 9600 > 9000 = (p²−1)p·E`)。 -/
theorem witness_five_exceeds_pair_budget :
    ¬ ((5:ℝ) ^ ((24:ℝ) / 75) ≤ ∏ i ∈ Finset.Icc 1 2, axDecay 5 i) := by
  haveI : Fact (Nat.Prime 5) := ⟨Nat.prime_five⟩
  intro h
  have h2 := (pairBudget_iff (p := 5) (J := 24) (e := 3) (by norm_num)).mp (by
    norm_num at h ⊢; exact h)
  norm_num [pairBudget] at h2

/-! ## §6 `.src`(原典の対応箇所) -/

def defect_eq.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def jump_add_defect_le.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def min_cost_fits_of_le_three.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def witness_gap_exceeds.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def pairBudget_iff.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end JumpDefect

end ABC3.Found.PGC

section Audit
open ABC3.Found.PGC
#print axioms JumpDefect.defect_eq
#print axioms JumpDefect.jump_add_defect
#print axioms JumpDefect.jump_add_defect_le
#print axioms JumpDefect.jump_add_defect_eq_iff
#print axioms JumpDefect.defect_nonneg_iff
#print axioms JumpDefect.contrCost_fits
#print axioms JumpDefect.projCost_eq_of_stable
#print axioms JumpDefect.contrCost_fits_of_stable
#print axioms JumpDefect.projCost_fits_iff_of_stable
#print axioms JumpDefect.projCost_fits_of_le_three
#print axioms JumpDefect.min_cost_fits_of_le_three
#print axioms JumpDefect.covered_iff
#print axioms JumpDefect.covered_of_le_three
#print axioms JumpDefect.gap_pos_of_five_le
#print axioms JumpDefect.witness_five_admissible
#print axioms JumpDefect.witness_five_costs
#print axioms JumpDefect.witness_five_exceeds_budget
#print axioms JumpDefect.witness_five_exceeds_ax_constant
#print axioms JumpDefect.witness_gap_admissible
#print axioms JumpDefect.witness_gap_projCost
#print axioms JumpDefect.witness_gap_exceeds
#print axioms JumpDefect.pairBudget_iff
#print axioms JumpDefect.Zeta27.costs
#print axioms JumpDefect.Zeta27.proj_fits
#print axioms JumpDefect.Zeta27.zeta81_proj_fits
#print axioms JumpDefect.slack_three_eq
#print axioms JumpDefect.witness_five_exceeds_pair_budget
end Audit
