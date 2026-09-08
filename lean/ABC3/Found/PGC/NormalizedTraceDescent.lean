import ABC3.Found.PGC.AxEpsilonDecay

/-!
# [pGC] 「正規化した跡」で `AxWildDescent K (axDecay p)` を出す道 —— ★勘定が合わない

持ち場が配ったのは「古典的な Ax–Sen–Tate は正規化した跡である。
`L/K` 有限次分離で `x ∈ L` のとき `y := Tr_{L/K}(x)/[L:K]` を近似元に取り、
`‖x − y‖` を different で押さえる」という道である。

★★**結論を先に書く: この道は `p ≥ 3` で閉じない。** 跡が与える指数は
sharp な指数のちょうど **`(p−1)` 倍**であり、`k = 1` の段ですでに
`axDecay p 1 = p^{1/(p−1)}` を超える。★`p = 2` に限れば `(p−1) = 1` なので
跡の上界は sharp と一致し、閉じる。★これが本ファイルの主結果 1 である。

★★**同じ勘定が「閉じる道」も特定した**(主結果 3、`FirstJumpRoute`)。
跡ではなく**第 1 跳び**で測ると総積がちょうど `p^{p/(p−1)²}` に収まる。

## §A 勘定(★Lean を書く前に手で閉じた。以下はその記録)

### A.1 跡が与える 1 段の損失

`L/M` を次数 `p` の完全分岐 Galois 層、下付き番号の跳びを `i`、
`e_L := e(L/ℚ_p)` とする。跡で近似元を作るには
`Tr_{L/M}(a) = 1` なる `a` を取り `x′ := Tr_{L/M}(a·x)` と置く。すると

  `x − Tr(a x) = Σ_{τ ∈ Gal(L/M)} τ(a)·(x − τ x)`   ⇒   `‖x − x′‖ ≤ ‖a‖·ε`

(★これが `TraceCore.norm_sub_traceAverage_le`。★分岐の語彙が 1 語も無い)。
`‖a‖` の最小値は Serre `Corps Locaux` V §3 Lemme 4
`Tr(𝔭_L^n) = 𝔭_M^{⌊(n+d)/e⌋}` から出る: `Tr(𝔭_L^n) = 𝒪_M` は
`−d ≤ n ≤ e−1−d` と同値なので最良は `v_L(a) = e − 1 − d`、すなわち

  `‖a‖ = p^{(d − e + 1)/e_L}`,  `d = (p−1)(i+1)`, `e = p`
  ⇒ ★`‖a‖ = p^{(p−1)·i / e_L}`.

### A.2 sharp な 1 段の損失(木がすでに持っている)

`RamificationJumpBound.norm_sub_digit_zero_le_rpow_mul` は
`d(x,M) = ‖π_L‖^{−i}·‖σx − x‖`、すなわち損失 `p^{i/e_L}` を与える。
★★**跡は指数を `(p−1)` 倍に膨らませる**(`JumpArith.traceLoss_eq_sharpLoss_rpow`)。

### A.3 `(p−1)i ≤ e_L` を入れて閉じた形にする

`RamificationJumpBound` の主結果 `(p−1)·i ≤ e_L` を入れると

* sharp: `p^{i/e_L} ≤ p^{1/(p−1)} = axDecay p 1`  ★ちょうど一致(等号が実現)
* 跡  : `p^{(p−1)i/e_L} ≤ p^{1} = p`         ★`k` に依らない定数

`p ≥ 3` では `axDecay p 1 = p^{1/(p−1)} < p` なので
★**跡の上界は `k = 1` ですでに目標を超える**(`JumpArith.axDecay_one_lt_natCast`)。
しかも定数 `p` を毎段払うと `∏_{k∈[1,n]} p = p^n` は非有界
(`JumpArith.constantLoss_prod_unbounded`)なので、
★**塔で積を取っても `p^{p/(p−1)²}` に収まらない。**

### A.4 反例データ `p = 3`, `ℚ₃(ζ₂₇)/ℚ₃(ζ₃)`(持ち場の指示どおり数値を入れた)

`L = ℚ₃(ζ₂₇)`, `F = ℚ₃(ζ₃)`, `[L:F] = 9`, `e_L = 18`, wild 深さ `k = 2`。
`Gal(L/ℚ₃)` の下付き分岐群は `G_u = Gal(L/ℚ₃(ζ_{3^j}))` (`3^{j−1} ≤ u ≤ 3^j − 1`) なので
`H := Gal(L/F)` は `H_u = H` (`u ≤ 2`)、`H_u = Gal(L/ℚ₃(ζ₉))` (`3 ≤ u ≤ 8`)、`H_9 = 1`。

| 測り方 | 使う跳び | 指数 | 値 | `axDecay 3 2 = 3^{1/6}` と比較 |
|---|---|---|---|---|
| 跡(上の層 `L/ℚ₃(ζ₉)`) | `i = 8` | `(3−1)·8/18 = 8/9` | `3^{0.889}` | ★超える(`axDecay_three_two_lt_traceBound`) |
| 跡(`L/F` 全体、`d_{L/F} = 36`) | —— | `(36−9+1)/18 = 14/9` | `3^{1.556}` | ★もっと悪い |
| sharp(上の層) | `i = 8` | `8/18 = 4/9` | `3^{0.444}` | ★これでも超える |
| ★★sharp(**第 1 跳び**) | `i₁ = 2` | `2/18 = 1/9` | `3^{0.111}` | ★★収まる(`sharp_three_two_le_axDecay`) |

`k = 1` の層 `ℚ₃(ζ₉)/ℚ₃(ζ₃)` (`e = 6`, 跳び `i = 2`, `d = (3−1)(2+1) = 6`) でも同じ:
跡は `3^{(6−3+1)/6} = 3^{2/3}`、sharp は `3^{2/6} = 3^{1/3}`、目標は `3^{1/2}`。
★跡は超え(`axDecay_three_one_lt_traceBound`)、sharp は収まる(`sharp_three_one_le_axDecay`)。

★★この `i = 2` は手で 2 度計算した。`σ(ζ₉) = ζ₉⁴`, `π = ζ₉ − 1` として
`σπ/π − 1 = 3 + 6π + 4π² + π³`。`v(3) = 6`, `v(4π²) = 2`, `v(π³) = 3` なので `v = 2`。
★1 度目は `4π²` を落として `v = 3` と誤り、`i = 3` を出した。
標準結果(`G_u = Gal(L/ℚ_p(ζ_{p^j}))` の区間)と突き合わせて発見した。

## §B 「閉じる道」—— 第 1 跳びで測る(★本ファイルの主結果 3)

A.4 の表の最終行が偶然でないことは Herbrand で説明できる:

`L/F` を完全分岐 Galois、`i₁` を `Gal(L/F)` の下付き第 1 跳び(`L` の番号)、
`E₁` を `G_{i₁+1}` の固定体(`[E₁:F] = p`)とすると、
`Gal(L/E₁) = G_{i₁+1}` の跳びはすべて `i₁ + 1` 以上なので
`φ_{L/E₁}(u) = u` (`u ≤ i₁+1`)、したがって `ψ_{L/E₁}(j) = j`。
Herbrand `(G/H)_{φ(u)} = G_u H/H` より

  ★`i₁ = j`  (`j` は**下の層** `E₁/F` の、`E₁` 自身の番号での跳び)

`RamificationJumpBound` を `E₁/F` に当てると `(p−1)·j ≤ e_{E₁}` であり、
`e_L = e(L/E₁)·e_{E₁}` かつ `p^{k−1} ∣ e(L/E₁)` なので

  `i₁/e_L = j/(e(L/E₁)·e_{E₁}) ≤ (e_{E₁}/(p−1))/(p^{k−1}·e_{E₁}) = (1/(p−1))·p^{1−k}`

  ⇒ ★★`p^{i₁/e_L} ≤ axDecay p k`  ★★**ちょうど一致する**。

この不等式が `JumpArith.rpow_div_le_axDecay` であり、
それを `AxWildDescent` に差し込んだのが `FirstJumpRoute.axLemma_of_firstJump` /
`axSenTate_of_firstJump` である。★**次の波が埋めるべきものは
「wild 深さ `k` の `x` に対し `‖x − x′‖ ≤ p^{j/(m·e)}·ε` を満たす `x′`(深さ `< k`)を作る」
ただ 1 点**で、`j`(下の層の跳び)・`m`(`p^{k−1} ≤ m`)・`e`(`(p−1)j ≤ e`)は
その 3 つの仮定を満たしさえすればよい。分岐の語彙は残っていない。

## §C 在庫の測定(★自分でコマンドを叩いた)

```
grep -n "norm_sum_le_of_forall_le" .cache/mathlib-index.txt      → ★0 件
  ★しかし #check @IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg は在る:
     0 ≤ C → (∀ i ∈ s, ‖f i‖ ≤ C) → ‖∑ i ∈ s, f i‖ ≤ C
  ★`AxEpsilonDecay.lean` の docstring が「索引の嘘」と書いたのと同じ穴。
grep -n "div_le_div_iff" .cache/mathlib-index.txt
  → `div_le_div_iff` は無い(`div_le_div_iff₀` / `_of_pos_left` / `_of_pos_right`)。
     素で書くと `Unknown identifier` になる。★`lean-idioms.md` 2116 行に既出。
grep -n "MulSemiringAction" .cache/mathlib-index.txt
  → `MulSemiringAction` は `DistribMulAction` の拡張。`smul_mul'` は
     `MulDistribMulAction` 用だが `MulSemiringAction` から解決した(推論が通した)。
grep -rn "intTrace|Algebra.trace" lean/ABC3/Found/PGC/ --include=*.lean
  → ★本ファイル自身の docstring 1 行だけ(＝本ファイル以前は 0 件)。
  ★持ち場の見立て(「PGC 側に跡写像の宣言が 0 件、未着手の道」)は正しかった。
  ★ただし未着手なのは「使えない」からではなく「使っても足りない」からだと本ファイルが測った。
  ★同じ grep を `lean/ABC3/Found/` で撃つと `Falt1/AlmostDerivation.lean` ほか 10 本以上
    出る(`GenEll/DifferentTameGlobal.lean` に `Algebra.intTrace`)。★在庫はあるが要らない。
```

## §D 逸脱の記録

* ★**具体層(`K.closure` への代入)を作らなかった。** 理由は §A.3:
  跡の抽象核は `sorry` 0 で証明したが、それを `PAdicLocalField` に代入しても
  `p ≥ 3` では `axDecay p 1` を超えるので `AxWildDescent K (axDecay p)` は出ない。
  ★代入に必要な配管(剰余類代表の `Finset (K.absGal)`、有限層の構成)は
  結論が出ないと分かっている道に払う費用なので払っていない。
* ★`p = 2` に限れば跡で閉じる(`traceExponent_eq_sharp_two` で `(p−1) = 1`)。
  ただし `AxSenTate` は全 `p` で要るので単独では役に立たない。
* ★`FirstJumpRoute` の 2 本は `AxWildDescent K c` を仮定として残しており、
  **`AxWildDescent` 自体は埋めていない**。埋めたのは「その `c` が `axDecay p` に
  収まるための閉じた十分条件」である。
-/

open Finset

namespace ABC3.Found.PGC

/-! ## §1 抽象核 —— 正規化した跡

★以下の 3 本には**分岐・付値・Galois・`p` 進の語彙が 1 語も出てこない**。
超距離の付いた可換環と、それに環自己同型として等長に作用する群だけである。 -/

namespace TraceCore

variable {R : Type*} [NormedCommRing R] [IsUltrametricDist R]

/-- ★★**抽象核 1(重み付き平均)** —— 重みの総和が `1` なら、近似の誤差は
「重みの上界 × 各点の誤差」で押さえられる。

`x − Σ a_i z_i = Σ a_i (x − z_i)` と超距離の一様上界だけ。
★★これが「正規化した跡」の全部であり、跡・分岐・付値は一切要らない。 -/
theorem norm_sub_sum_mul_le {ι : Type*} {s : Finset ι} {a z : ι → R} {x : R} {M ε : ℝ}
    (hM : 0 ≤ M) (hε : 0 ≤ ε) (hsum : ∑ i ∈ s, a i = 1)
    (ha : ∀ i ∈ s, ‖a i‖ ≤ M) (hz : ∀ i ∈ s, ‖x - z i‖ ≤ ε) :
    ‖x - ∑ i ∈ s, a i * z i‖ ≤ M * ε := by
  have key : x - ∑ i ∈ s, a i * z i = ∑ i ∈ s, a i * (x - z i) := by
    simp only [mul_sub]
    rw [Finset.sum_sub_distrib, ← Finset.sum_mul, hsum, one_mul]
  rw [key]
  refine IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg (by positivity) ?_
  intro i hi
  exact le_trans (norm_mul_le _ _) (mul_le_mul (ha i hi) (hz i hi) (norm_nonneg _) hM)

def norm_sub_sum_mul_le.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★**抽象核 2(正規化した跡)** —— 群 `G` が環 `R` に等長な環自己同型として作用し、
`Σ_{g ∈ s} g • a = 1`(★これが「跡が `1` になる元」)であるとき

  `‖x − Σ_{g ∈ s} g • (a·x)‖ ≤ ‖a‖ · max_g ‖x − g • x‖`.

★古典的な Ax–Sen–Tate の「`y := Tr_{L/K}(a x)`」がこれである。
★`a = 1/n`(`n = |G|`)を入れると重心(`norm_sub_barycenter_le`)に退化する。
★★**分岐・付値・Galois の語彙が 1 語も無い。** -/
theorem norm_sub_traceAverage_le {G : Type*} [Group G] [MulSemiringAction G R]
    (hiso : ∀ (g : G) (r : R), ‖g • r‖ = ‖r‖)
    {s : Finset G} {a x : R} {ε : ℝ} (hε : 0 ≤ ε)
    (hsum : ∑ g ∈ s, g • a = 1) (hx : ∀ g ∈ s, ‖x - g • x‖ ≤ ε) :
    ‖x - ∑ g ∈ s, g • (a * x)‖ ≤ ‖a‖ * ε := by
  have hrw : ∀ g ∈ s, g • (a * x) = (g • a) * (g • x) := fun g _ => smul_mul' g a x
  rw [Finset.sum_congr rfl hrw]
  exact norm_sub_sum_mul_le (norm_nonneg a) hε hsum (fun g _ => le_of_eq (hiso g a)) hx

def norm_sub_traceAverage_le.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- **抽象核 3(重心)** —— `(|s| : R)·a = 1` の場合。
★木の `exists_norm_sub_algebraMap_le_div_norm_natDegree`(定数 `|n|^{−1}`)の抽象形。 -/
theorem norm_sub_barycenter_le {ι : Type*} {s : Finset ι} {z : ι → R} {x a : R} {ε : ℝ}
    (hε : 0 ≤ ε) (hna : (s.card : R) * a = 1) (hz : ∀ i ∈ s, ‖x - z i‖ ≤ ε) :
    ‖x - ∑ i ∈ s, a * z i‖ ≤ ‖a‖ * ε := by
  refine norm_sub_sum_mul_le (norm_nonneg a) hε ?_ (fun _ _ => le_refl _) hz
  rw [Finset.sum_const, nsmul_eq_mul]
  exact hna

def norm_sub_barycenter_le.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end TraceCore

/-! ## §2 指数の勘定 —— 抽象核(ℕ と ℝ だけ)

★ここにも分岐・付値・Galois は出てこない。出るのは `axDecay` の定義だけである。 -/

namespace JumpArith

variable {p : ℕ} [Fact p.Prime]

/-- ★★★**閉じる十分条件**(本ファイルの主結果 3、§B の勘定)。

`p^{j/(m·e)} ≤ axDecay p k`  を、次の 3 つだけから出す:

* `p^{k−1} ≤ m`  (★塔の下から `k−1` 段ぶんの分岐指数が `m` に入っている)
* `0 < e`
* `(p−1)·j ≤ e`  (★`RamificationJumpBound` の主結果を**下の層**に当てた形)

★★分岐の言葉に戻すと `j` = 下の層 `E₁/F` の跳び、`m` = `e(L/E₁)`、`e` = `e_{E₁}`、
`m·e = e_L` である。★★これが `AxWildDescent K (axDecay p)` に残っている
ただ 1 点を「跡」ではなく「第 1 跳び」で測る道の勘定の全部。 -/
theorem rpow_div_le_axDecay {k j m e : ℕ}
    (hm : p ^ (k - 1) ≤ m) (he : 0 < e) (hjump : (p - 1) * j ≤ e) :
    (p : ℝ) ^ ((j : ℝ) / ((m : ℝ) * (e : ℝ))) ≤ axDecay p k := by
  have hp2 : 2 ≤ p := (Fact.out : p.Prime).two_le
  have h1 : (1:ℝ) < (p : ℝ) := by exact_mod_cast (Fact.out : p.Prime).one_lt
  have hjumpR : ((p : ℝ) - 1) * (j : ℝ) ≤ (e : ℝ) := by
    have := (Nat.cast_le (α := ℝ)).mpr hjump
    rwa [Nat.cast_mul, Nat.cast_sub (by omega), Nat.cast_one] at this
  have hmR : (p : ℝ) ^ (k - 1) ≤ (m : ℝ) := by
    have := (Nat.cast_le (α := ℝ)).mpr hm
    rwa [Nat.cast_pow] at this
  have hpk : (0:ℝ) < (p : ℝ) ^ (k - 1) := by positivity
  have hm0 : (0:ℝ) < (m : ℝ) := lt_of_lt_of_le hpk hmR
  have he0 : (0:ℝ) < (e : ℝ) := by exact_mod_cast he
  rw [axDecay]
  refine Real.rpow_le_rpow_of_exponent_le (le_of_lt h1) ?_
  have hrhs : (1 / ((p:ℝ) - 1)) * (1 / (p:ℝ)) ^ (k - 1)
      = 1 / (((p:ℝ) - 1) * (p:ℝ) ^ (k - 1)) := by
    rw [one_div_pow, div_mul_div_comm, one_mul]
  rw [hrhs, div_le_div_iff₀ (by positivity) (by positivity)]
  have step : (j:ℝ) * (((p:ℝ) - 1) * (p:ℝ) ^ (k-1))
      = (((p:ℝ) - 1) * (j:ℝ)) * (p:ℝ) ^ (k-1) := by ring
  rw [step, one_mul]
  calc (((p:ℝ) - 1) * (j:ℝ)) * (p:ℝ) ^ (k-1) ≤ (e:ℝ) * (p:ℝ) ^ (k-1) :=
        mul_le_mul_of_nonneg_right hjumpR (le_of_lt hpk)
    _ ≤ (e:ℝ) * (m:ℝ) := mul_le_mul_of_nonneg_left hmR (le_of_lt he0)
    _ = (m:ℝ) * (e:ℝ) := mul_comm _ _

def rpow_div_le_axDecay.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★**主結果 1 —— 跡は sharp の `(P−1)` 乗である。**

different が与える損失 `P^{(P−1)i/e}` は、`RamificationJumpBound` が与える
sharp な損失 `P^{i/e}` の **`(P−1)` 乗**にちょうど等しい。

★★`P = 2` なら両者は一致する(だから `p = 2` では跡の道も閉じる)。
★★`P ≥ 3` では跡が真に大きい —— これが本ファイルが「閉じない」と測った理由である。 -/
theorem traceLoss_eq_sharpLoss_rpow (P i e : ℝ) (hP : 0 ≤ P) :
    P ^ (((P - 1) * i) / e) = (P ^ (i / e)) ^ (P - 1) := by
  rw [← Real.rpow_mul hP]
  congr 1
  ring

def traceLoss_eq_sharpLoss_rpow.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- `p = 2` なら跡の指数と sharp な指数は一致する。 -/
theorem traceExponent_eq_sharp_two (i e : ℝ) : (((2:ℝ) - 1) * i) / e = i / e := by ring

def traceExponent_eq_sharp_two.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★**主結果 1 の要**: `p ≥ 3` では `axDecay p 1 < p`。

跡の道が `(p−1)i ≤ e_L` から出せる上界は `p^{(p−1)i/e_L} ≤ p` であり、
★これは `k = 1` の段ですでに目標 `axDecay p 1 = p^{1/(p−1)}` を超える。 -/
theorem axDecay_one_lt_natCast (hp : 3 ≤ p) : axDecay p 1 < (p : ℝ) := by
  have h1 : (1:ℝ) < (p : ℝ) := by exact_mod_cast (Fact.out : p.Prime).one_lt
  have h3 : (3:ℝ) ≤ (p : ℝ) := by exact_mod_cast hp
  rw [axDecay_one]
  have hlt : (p:ℝ) ^ (1 / ((p:ℝ) - 1)) < (p:ℝ) ^ (1:ℝ) := by
    refine Real.rpow_lt_rpow_of_exponent_lt h1 ?_
    rw [div_lt_one (by linarith)]
    linarith
  rwa [Real.rpow_one] at hlt

def axDecay_one_lt_natCast.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★**定数の損失は塔で発散する** —— 毎段 `p` を払うと `∏_{[1,n]} p = p^n` は非有界。

跡の道は(`(p−1)i ≤ e_L` しか使えないので)`k` に依らない定数 `p` しか出せない。
★★したがって塔で積を取っても `axConstant p = p^{p/(p−1)²}` に収まらない。 -/
theorem constantLoss_prod_unbounded (C : ℝ) :
    ∃ n : ℕ, C < ∏ _k ∈ Finset.Icc 1 n, (p:ℝ) := by
  have h1 : (1:ℝ) < (p : ℝ) := by exact_mod_cast (Fact.out : p.Prime).one_lt
  exact prod_unbounded_of_one_lt h1 (c := fun _ => (p:ℝ)) (fun _ => le_refl _) C

def constantLoss_prod_unbounded.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ### §2.1 反例データ `p = 3` の数値(§A.4 の表を機械に検査させたもの) -/

theorem axDecay_three_one_eq : axDecay 3 1 = (3:ℝ) ^ ((1:ℝ)/2) := by
  norm_num [axDecay]

def axDecay_three_one_eq.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

theorem axDecay_three_two_eq : axDecay 3 2 = (3:ℝ) ^ ((1:ℝ)/6) := by
  norm_num [axDecay]

def axDecay_three_two_eq.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★`ℚ₃(ζ₉)/ℚ₃(ζ₃)`(`e = 6`, 跳び `i = 2`, different 指数 `d = 6`)で
跡が出す上界 `3^{(6−3+1)/6} = 3^{2/3}` は目標 `axDecay 3 1 = 3^{1/2}` を**超える**。 -/
theorem axDecay_three_one_lt_traceBound : axDecay 3 1 < (3:ℝ) ^ ((2:ℝ)/3) := by
  rw [axDecay_three_one_eq]
  exact (Real.rpow_lt_rpow_left_iff (by norm_num)).mpr (by norm_num)

def axDecay_three_one_lt_traceBound.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- 同じ層の sharp な損失 `3^{2/6} = 3^{1/3}` は目標に**収まる**。
★跡と sharp の差(指数で `2` 倍 `= p−1`)がそのまま成否を分けている。 -/
theorem sharp_three_one_le_axDecay : (3:ℝ) ^ ((1:ℝ)/3) ≤ axDecay 3 1 := by
  rw [axDecay_three_one_eq]
  exact (Real.rpow_le_rpow_left_iff (by norm_num)).mpr (by norm_num)

def sharp_three_one_le_axDecay.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★持ち場が名指しした反例データ `p = 3`, `ℚ₃(ζ₂₇)/ℚ₃(ζ₃)`(`i_m = 8`, `e_L = 18`)。
上の層に跡を当てると `3^{(3−1)·8/18} = 3^{8/9}`、目標は `axDecay 3 2 = 3^{1/6}`。
★★**指数で 5.33 倍超える。** -/
theorem axDecay_three_two_lt_traceBound : axDecay 3 2 < (3:ℝ) ^ ((8:ℝ)/9) := by
  rw [axDecay_three_two_eq]
  exact (Real.rpow_lt_rpow_left_iff (by norm_num)).mpr (by norm_num)

def axDecay_three_two_lt_traceBound.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★同じデータを**第 1 跳び**で測ると収まる: `j = 2`(下の層 `ℚ₃(ζ₉)/ℚ₃(ζ₃)` の跳び)、
`m = 3 = e(L/ℚ₃(ζ₉)) = p^{k−1}`、`e = 6 = e_{ℚ₃(ζ₉)}`、`m·e = 18 = e_L`。
`(p−1)j = 4 ≤ 6 = e` なので `rpow_div_le_axDecay` が通り
`3^{2/18} = 3^{1/9} ≤ 3^{1/6} = axDecay 3 2`。 -/
theorem sharp_three_two_le_axDecay : (3:ℝ) ^ ((1:ℝ)/9) ≤ axDecay 3 2 := by
  haveI : Fact (Nat.Prime 3) := ⟨Nat.prime_three⟩
  have h := rpow_div_le_axDecay (p := 3) (k := 2) (j := 2) (m := 3) (e := 6)
    (by norm_num) (by norm_num) (by norm_num)
  have hexp : ((2:ℕ) : ℝ) / (((3:ℕ) : ℝ) * ((6:ℕ) : ℝ)) = (1:ℝ)/9 := by norm_num
  rw [hexp] at h
  simpa using h

def sharp_three_two_le_axDecay.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end JumpArith

/-! ## §3 具体層 —— 第 1 跳びの道を `AxWildDescent` に差し込む -/

namespace FirstJumpRoute

open ABC3.Skeleton.PGC

variable {p : ℕ} [Fact p.Prime]

/-- ★★★**主結果 3** —— 1 段の損失が「第 1 跳び `j` / (`m·e`)」の形で押さえられ、
`p^{k−1} ≤ m` と `(p−1)j ≤ e` が成り立つなら **`AxLemma K (axConstant p)`**。

★★これが `AxWildDescent K (axDecay p)` を出す**閉じた十分条件**である。
残っているのは「そういう `x′` を作る」具体層だけで、
★**定数の勘定はここで完全に閉じている**(`JumpArith.rpow_div_le_axDecay`)。 -/
theorem axLemma_of_firstJump (K : PAdicLocalField p) {c : ℕ → ℝ} {j m e : ℕ → ℕ}
    (hc : ∀ k, 1 ≤ c k)
    (hform : ∀ k, c k ≤ (p:ℝ) ^ ((j k : ℝ) / ((m k : ℝ) * (e k : ℝ))))
    (hm : ∀ k, p ^ (k - 1) ≤ m k) (he : ∀ k, 0 < e k)
    (hjump : ∀ k, (p - 1) * j k ≤ e k)
    (h : AxWildDescent K c) : AxLemma K (axConstant p) :=
  axLemma_of_axDecay K hc
    (fun k => (hform k).trans (JumpArith.rpow_div_le_axDecay (hm k) (he k) (hjump k))) h

def axLemma_of_firstJump.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★そこから `AxSenTate K`。 -/
theorem axSenTate_of_firstJump (K : PAdicLocalField p) {c : ℕ → ℝ} {j m e : ℕ → ℕ}
    (hc : ∀ k, 1 ≤ c k)
    (hform : ∀ k, c k ≤ (p:ℝ) ^ ((j k : ℝ) / ((m k : ℝ) * (e k : ℝ))))
    (hm : ∀ k, p ^ (k - 1) ≤ m k) (he : ∀ k, 0 < e k)
    (hjump : ∀ k, (p - 1) * j k ≤ e k)
    (h : AxWildDescent K c) : AxSenTate K :=
  axSenTate_of_axLemma (le_trans zero_le_one (one_le_axConstant p))
    (axLemma_of_firstJump K hc hform hm he hjump h)

def axSenTate_of_firstJump.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end FirstJumpRoute

end ABC3.Found.PGC

section Audit
open ABC3.Found.PGC
#print axioms TraceCore.norm_sub_sum_mul_le
#print axioms TraceCore.norm_sub_traceAverage_le
#print axioms TraceCore.norm_sub_barycenter_le
#print axioms JumpArith.rpow_div_le_axDecay
#print axioms JumpArith.traceLoss_eq_sharpLoss_rpow
#print axioms JumpArith.traceExponent_eq_sharp_two
#print axioms JumpArith.axDecay_one_lt_natCast
#print axioms JumpArith.constantLoss_prod_unbounded
#print axioms JumpArith.axDecay_three_one_lt_traceBound
#print axioms JumpArith.axDecay_three_two_lt_traceBound
#print axioms JumpArith.sharp_three_one_le_axDecay
#print axioms JumpArith.sharp_three_two_le_axDecay
#print axioms FirstJumpRoute.axLemma_of_firstJump
#print axioms FirstJumpRoute.axSenTate_of_firstJump
end Audit
