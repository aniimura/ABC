import ABC3.Meta.Claim
import ABC3.Interface.GenEll.EllModuli
import ABC3.Found.GenEll.Elementary
import ABC3.Found.GenEll.PrimesOfSize
import ABC3.Found.GenEll.PrimeNumberTheorem
import ABC3.Found.GenEll.PrimeConstants
import ABC3.Skeleton.GenEll.GaloisImage
import Mathlib.NumberTheory.Chebyshev

/-!
# [GenEll] §4 Primes of Prescribed Size —— 5 件(`Skeleton`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、
物理 p.20–p.23。**260 dpi 目視確認 2026-08-16**。

原文 (GenEll p.20):
> Lemma 4.1. (The Existence of Primes of Prescribed Size) Write

`Corollary 4.4` が**本トラックの北極星**である(`ResearchPaper/genell-goal.md`)。

## ★★実測: `Chebyshev.theta` は原文の `θ` と `ℝ′_{>0}` 上で**一致する**

mathlib の `Chebyshev.theta x = Σ_{p ≤ x, p 素数} log p`(`⌊x⌋₊` を経由)。
原文の `θ(x) = Σ_{p < x} log p` は **`x ∈ ℝ′_{>0}`(素数でない正実数)** に限って使われる。

★`x` が素数でない正実数なら `p ≤ x` と `p < x` は**同じ素数の集合を与える**
(`p = x` は `x` が素数のときしか起きない)。
★★**したがって原文が `ℝ′_{>0}` を導入したのは偶然ではなく、
この一致を成り立たせるためだと読める。** よって mathlib の `Chebyshev.theta` を
そのまま使ってよい——posit する必要がない。

## ★★実測: 素数定理は mathlib に**無い**

`Mathlib/NumberTheory/Chebyshev.lean` は Chebyshev 型の評価
(`theta_ge'` / `psi_ge'` / `pi_le_log4_mul_div` など)を持つが、
**`θ(x)/x → 1`(素数定理)は無い**(2026-08-16、同ファイルの全定理名を確認)。

★これは `Remark 4.1.1` の後半にそのまま効く——
★しかし **`Lemma 4.1` 自身には効かない**。`Lemma 4.1` は条件 (i)(ii) を
**仮説として受け取る**形になっているからである。
`Remark 4.1.1` が「条件 (i) は素数定理の帰結」と書いているのは、
**素数定理を圏外に置く**ためのただし書きである。

## ★2 理論なしで届くもの

| 項目 | Arakelov | Galois 表現 | 実際に要るもの |
|---|---|---|---|
| `Lemma 4.1` | 要らない | 要らない | ★**背理法と 6 行の不等式計算だけ** |
| `Lemma 4.2` | 要らない | 要らない | ★**`log(H+1) ≤ (3H/2)log 2`(初等)** |
| `Remark 4.1.1` | 要らない | 要らない | ★**素数定理**(mathlib に無い) |
| `Corollary 4.3` | 要る | 要る | `Theorem 3.8` |
| `Corollary 4.4` | 要る | 要る | `Theorem 3.8` |

## ★★★★★★★★ 2026-08-26——§4 の `sorry` が 0 になった(第 367 ブロック)

`Corollary 4.3` / `Corollary 4.4` を `Theorem 3.8` から**導出**した。

### ★★★実際に使ったもの(すべて原文 p.22-23 が引いている)

| 段 | どこから |
|---|---|
| `l` が不分岐 ⇒ `SL₂ ⊆ 像` と全射は同じ | 界面 `imageSurjective_of_containsSL2`(原文 p.22 の冒頭) |
| `x_{S∘} ≤ x_S + (1 + 3/2)·23040d·deg∞` | **`Lemma 4.2`**(実装済み)+ 界面 `sum_localHt_eq` |
| `deg∞ ≤ 12·(6/5)·ht^Falt + C` | **`Proposition 3.4`**(第 361 で導出済み) |
| `x_ϵ ≤ x_S` にしてよい | **`exists_finset_primes_sum_log_gt`**(第 366) |
| 条件 (i)(ii) を満たす定数 | **`exists_cond_i_ii`**(第 366、`ϵ = 1/6`) |
| `M = 1` の素数の存在 | **`Lemma 4.1`**(実装済み) |

### ★★★★★係数 900 と 100 の差は `h` だけである

原文 p.23 の帳尻は `2·3·12 + 8·100 ≤ 100 + 800 = 900`:

    Corollary 4.3: h = 23040·100d·(ht^Falt + C′·d^ϵ)  ⇒  8h が 800 を出す  ⇒ 900、末項は C·d^{1+ϵ}
    Corollary 4.4: h = 0                                ⇒  8h は消える    ⇒ 100、末項は C·d

★どちらも `2·x_bad` から来る `2·3·12 = 72` が土台である。
★★`ht^Falt < 0` のとき `72 ≤ 100` は逆向きになるが、
差 `28·23040·d·(ht^Falt + B) ≥ 0` がちょうど埋める。

### ★★★★★★`30` と素という条件の扱い方

`Corollary 4.4` は `Theorem 3.8` の**条件 (b)** を使うので `l` が `30` と素でなければならない。
★**除外集合に `{2, 3, 5}` を入れる**だけですむ——その分の `log 30` は定数なので
`C·d` に吸収される(`d ≥ 1`)。

### 逸脱

`Corollary 4.3` の仮説 `MinimalField E` は**使っていない**(弱めてはいない)。
★原文は `d = [L:ℚ]` を最小定義体の次数として使うが、上の導出は
どの定義体でも通る——仮説を残したのは statement を原文に合わせるためである。
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta ABC3.Interface.GenEll

/-! ## Lemma 4.1 —— 指定された大きさの素数の存在 -/

/-- **[GenEll] Lemma 4.1**(The Existence of Primes of Prescribed Size)。

原文 (GenEll p.20):
> Lemma 4.1. (The Existence of Primes of Prescribed Size) Write

`ℝ′_{>0}` を `ℝ_{>0}` から素数を除いたもの、`θ(x) ≝ Σ_{p<x} log(p)` とする。
`M` を正整数、`ϵ, x_ϵ, C_ϵ ∈ ℝ_{>0}` を `0 < ϵ < 1/4`、`ϵ·x_ϵ > C_ϵ` かつ
(i) `(5/4)·x + C_ϵ > θ(x)`(全ての `x ∈ ℝ′_{>0}`)、
`θ(x) > (1−ϵ)x`(`x ≥ x_ϵ` なる `x ∈ ℝ′_{>0}`)、
(ii) `M·log(x) ≤ ϵ·x`(`x ≥ x_ϵ` なる `x ∈ ℝ_{>0}`)
を満たすものとする。`x_A ≝ Σ_{p∈A} log(p)` と書く。

このとき**任意の非負 `h ∈ ℝ` と `x_A > x_ϵ` なる任意の有限素数集合 `A` に対し、
`M` 個の相異なる素数 `p_1,…,p_M` が存在して `p_j ∉ A` かつ
`h ≤ p_j ≤ (1+6ϵ)·x_A + 8h`**。

★★**条件 (i)(ii) は仮説である。** 補題は「それらを満たす `ϵ, x_ϵ, C_ϵ` が
**与えられたとき**」を主張しており、**存在することは主張していない**。
★ゆえに素数定理を証明する必要がない——`Remark 4.1.1` がそれを明言する。

★`Chebyshev.theta` をそのまま使えることは上の docstring で実測した。

★★**本 statement は `sorry` ではない**——`Found/GenEll/PrimesOfSize.lean` の
実装をそのまま参照している。★原文が末尾で「WLOG `y_A, (1+δ)h ∈ ℝ′_{>0}`」と
済ませた段は、**`θ` の側を `ℝ_{>0}` 全体へ延ばす**ことで消してある
(`theta_le_of_cond_i` / `theta_gt_of_cond_ii`)。 -/
theorem lemma_4_1 (M : ℕ) (hM : 0 < M)
    (eps xeps Ceps : ℝ) (heps : 0 < eps) (hxeps : 0 < xeps) (hCeps : 0 < Ceps)
    (heps4 : eps < 1 / 4) (hxC : Ceps < eps * xeps)
    (hi1 : ∀ x : ℝ, 0 < x → (∀ p : ℕ, p.Prime → (p : ℝ) ≠ x) →
      Chebyshev.theta x < 5 / 4 * x + Ceps)
    (hi2 : ∀ x : ℝ, 0 < x → (∀ p : ℕ, p.Prime → (p : ℝ) ≠ x) → xeps ≤ x →
      (1 - eps) * x < Chebyshev.theta x)
    (hii : ∀ x : ℝ, 0 < x → xeps ≤ x → (M : ℝ) * Real.log x ≤ eps * x)
    (h : ℝ) (hh : 0 ≤ h) (A : Finset ℕ) (hA : ∀ p ∈ A, p.Prime)
    (hAx : xeps < ∑ p ∈ A, Real.log p) :
    ∃ P : Finset ℕ, P.card = M ∧ (∀ p ∈ P, p.Prime) ∧ (∀ p ∈ P, p ∉ A) ∧
      ∀ p ∈ P, h ≤ (p : ℝ) ∧
        (p : ℝ) ≤ (1 + 6 * eps) * (∑ q ∈ A, Real.log q) + 8 * h :=
  ABC3.Found.GenEll.lemma_4_1 M hM eps xeps Ceps heps hxeps hCeps heps4 hxC
    hi1 hi2 hii h hh A hA hAx

/-! ## Remark 4.1.1 —— 条件 (ii) は初等的、条件 (i) は素数定理 -/

/-- **[GenEll] Remark 4.1.1**。

原文 (GenEll p.21):
> which satisfy condition (ii) of Lemma 4.1 is entirely elementary. On the other hand, with regard to condition (i), the fact that

`Lemma 4.1` の条件 (ii) を満たす `ϵ, x_ϵ, C_ϵ` を見つけるのは**完全に初等的**。
一方、条件 (i) については `θ(x)/x → 1`(`x → ∞`)が
**素数定理のよく知られた帰結**である(cf. [Edw], p. 76)。

★★**この Remark が「素数定理は圏外」と明示している。**

★★**本 statement は `sorry` ではない**(2026-08-17 に閉じた)——
`Found/GenEll/PrimeNumberTheorem.lean` の `remark_4_1_1` を参照する。

- 第 1 項(条件 (ii) の初等性)は**我々が証明した**。
  原文が "entirely elementary" と書いた中身は `Real.isLittleO_log_id_atTop`
  (`log` は恒等写像より真に小さい)である。
- ★★第 2 項(`θ(x)/x → 1`)は**素数定理そのもの**で、**外部から借りた**——
  `PrimeNumberTheoremAnd.chebyshev_asymptotic`(`#print axioms` で
  標準 3 公理のみを実測。同 repo の `Wiener.lean` の sorry 2 件は経路外)。
  ★**原文も [Edw], p.76 を指すだけで証明していない。**
  借りたことは `Found` 側の docstring と `.src` に明記してある。 -/
theorem remark_4_1_1 :
    (∀ (M : ℕ) (eps : ℝ), 0 < eps →
        ∃ xeps : ℝ, 0 < xeps ∧ ∀ x : ℝ, xeps ≤ x → (M : ℝ) * Real.log x ≤ eps * x)
  ∧ Filter.Tendsto (fun x : ℝ => Chebyshev.theta x / x) Filter.atTop (nhds 1) :=
  ABC3.Found.GenEll.remark_4_1_1

/-! ## Lemma 4.2 —— いくつかの初等的な評価 -/

/-- **[GenEll] Lemma 4.2**(Some Elementary Estimates)。

原文 (GenEll p.21):
> Lemma 4.2. (Some Elementary Estimates) Let n be a positive integer;

`n` 正整数、`p_1,…,p_n` 素数、`h_1,…,h_n` 正整数、`h ≝ Σ h_j·log(p_j)` とすると
**`Σ log(p_j) ≤ h`** かつ **`Σ log(h_j) ≤ Σ log(h_j+1) ≤ 3h/2`**。

★★**純粋に初等的であり、今すぐ実装できる。**
原文の証明は 2 行——第 1・第 2 不等式は直ちに従い、第 3 不等式は
「任意の正整数 `H` について `log(H+1) ≤ (3H/2)·log(2)`」から従う(p.21 目視確認)。

★筋は完全に追える: `p_j ≥ 2` より `log(p_j) ≥ log 2` なので `h ≥ log(2)·Σ h_j`、
したがって `Σ log(h_j+1) ≤ (3/2)·log(2)·Σ h_j ≤ (3/2)·h`。

★★**本 statement は `sorry` ではない**——`Found/GenEll/Elementary.lean` の
実装をそのまま参照している。§4 の 5 項目のうち**実装まで済んでいる唯一の項目**である。 -/
theorem lemma_4_2 (n : ℕ) (hn : 0 < n) (p : Fin n → ℕ) (hp : ∀ j, (p j).Prime)
    (hh : Fin n → ℕ) (hhpos : ∀ j, 0 < hh j) :
    (∑ j, Real.log (p j)) ≤ (∑ j, (hh j : ℝ) * Real.log (p j))
  ∧ (∑ j, Real.log (hh j)) ≤ (∑ j, Real.log ((hh j : ℝ) + 1))
  ∧ (∑ j, Real.log ((hh j : ℝ) + 1))
      ≤ 3 / 2 * (∑ j, (hh j : ℝ) * Real.log (p j)) :=
  ABC3.Found.GenEll.lemma_4_2 n hn p hp hh hhpos

/-! ## Corollary 4.3 —— 退化する楕円曲線の完全 Galois 作用 -/

/-- **[GenEll] Corollary 4.3**(Full Galois Actions for Degenerating Elliptic Curves)。

原文 (GenEll p.22):
> Corollary 4.3. (Full Galois Actions for Degenerating Elliptic Curves)

`ϵ > 0` に対し**定数 `C > 0` と Galois-finite な `Exc` が存在して**次を満たす:
`E_L` を数体 `L` 上の楕円曲線で `L` が `[E_L]` の**最小定義体**、`[E_L] ∉ Exc`、
`S` を有限素数集合、`E_L` は**潜在的乗法還元の素点を 1 つ以上持つ**とすると、
`S` に属さない素数 `l_∘, l_•` が存在して:
(a) `l_∘, l_•` は潜在的乗法還元の素点および局所高さと素。さらに `l_•` は
`L` で分岐する `ℚ` の素数および分岐指数と素。
(b) `Gal(ℚ̄/L) → GL₂(ℤ_{l_∘})` の像が `SL₂(ℤ_{l_∘})` を含み、
`Gal(ℚ̄/L) → GL₂(ℤ_{l_•})` は**全射**。
(c) `l_∘ ≤ 23040·900d·ht^Falt + 2x_S + C·d^{1+ϵ}`、
`l_• ≤ 23040·900d·ht^Falt + 6d·log-diff + 2x_S + C·d^{1+ϵ}`。

★係数 `900` の由来は原文 p.23 に書かれている(`2·3·12 + 8·100 ≤ 100 + 800 = 900`)。 -/
theorem cor_4_3 (D : EllModuliData) (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, 0 < C ∧ ∃ Exc : Set D.EllClass, D.GaloisFinite Exc ∧
      ∀ (E : D.Curve) (S : Finset ℕ), (∀ p ∈ S, p.Prime) →
        D.MinimalField E → D.cls E ∉ Exc → D.HasPotMultRed E →
        ∃ lo lb : ℕ, Nat.Prime lo ∧ Nat.Prime lb ∧ lo ∉ S ∧ lb ∉ S
          ∧ (D.PrimeToMultPrimes E lo ∧ D.PrimeToLocalHeights E lo
              ∧ D.PrimeToMultPrimes E lb ∧ D.PrimeToLocalHeights E lb
              ∧ D.PrimeToRamification E lb)
          ∧ (D.ImageContainsSL2 E lo ∧ D.ImageSurjective E lb)
          ∧ (lo : ℝ) ≤ 23040 * 900 * (D.degOfDefinition E : ℝ) * D.faltingsHeight (D.cls E)
              + 2 * (∑ p ∈ S, Real.log p)
              + C * (D.degOfDefinition E : ℝ) ^ (1 + eps)
          ∧ (lb : ℝ) ≤ 23040 * 900 * (D.degOfDefinition E : ℝ) * D.faltingsHeight (D.cls E)
              + 6 * (D.degOfDefinition E : ℝ) * D.logDiffMell (D.cls E)
              + 2 * (∑ p ∈ S, Real.log p)
              + C * (D.degOfDefinition E : ℝ) ^ (1 + eps) := by
  classical
  obtain ⟨C₈, hC₈pos, Exc, hExc, h38⟩ := theorem_3_8 D ∅ D.compactlyBounded_empty eps heps
  obtain ⟨hdegbd, hhtbd, -, -⟩ := prop_3_4 D (1/5 : ℝ) (by norm_num)
  obtain ⟨A₁, hA₁⟩ := hdegbd
  obtain ⟨A₂, hA₂⟩ := hhtbd
  obtain ⟨Bl, hBl⟩ := D.faltingsHeight_bddBelow
  obtain ⟨xeps, Ceps, hxepos, hCepos, hCx, hi1, hi2, hii⟩ :=
    ABC3.Found.GenEll.exists_cond_i_ii 1
  obtain ⟨T, hTprime, hTgt⟩ := ABC3.Found.GenEll.exists_finset_primes_sum_log_gt xeps
  have hxTnn : (0:ℝ) ≤ ∑ p ∈ T, Real.log p :=
    Finset.sum_nonneg fun p _ => Real.log_natCast_nonneg p
  have hunion : ∀ U V : Finset ℕ,
      (∑ p ∈ U ∪ V, Real.log p) ≤ (∑ p ∈ U, Real.log p) + ∑ p ∈ V, Real.log p := by
    intro U V
    have hui : (∑ p ∈ U ∪ V, Real.log p) + ∑ p ∈ U ∩ V, Real.log p
        = (∑ p ∈ U, Real.log p) + ∑ p ∈ V, Real.log p := Finset.sum_union_inter
    have hnn : (0:ℝ) ≤ ∑ p ∈ U ∩ V, Real.log p :=
      Finset.sum_nonneg fun p _ => Real.log_natCast_nonneg p
    linarith
  refine ⟨23040 * (5 * |A₁ + A₂| + 800 * C₈ + 828 * |Bl|) + 2 * (∑ p ∈ T, Real.log p) + 1,
    by nlinarith [abs_nonneg (A₁ + A₂), abs_nonneg Bl, hC₈pos, hxTnn], Exc, hExc, ?_⟩
  intro E S hSprime _hmin hExcE hpot
  have hd1 : (1:ℝ) ≤ (D.degOfDefinition E : ℝ) := by exact_mod_cast D.degOfDefinition_pos E
  have hdpos : (0:ℝ) < (D.degOfDefinition E : ℝ) := by linarith
  have hdeps1 : (1:ℝ) ≤ (D.degOfDefinition E : ℝ) ^ eps := Real.one_le_rpow hd1 heps.le
  have hdd : (D.degOfDefinition E : ℝ) * (D.degOfDefinition E : ℝ) ^ eps
      = (D.degOfDefinition E : ℝ) ^ (1 + eps) := by
    rw [Real.rpow_add hdpos, Real.rpow_one]
  have hdP : (D.degOfDefinition E : ℝ) ≤ (D.degOfDefinition E : ℝ) ^ (1 + eps) := by
    rw [← hdd]; nlinarith
  have hP1 : (1:ℝ) ≤ (D.degOfDefinition E : ℝ) ^ (1 + eps) := le_trans hd1 hdP
  have hBnn : (0:ℝ) ≤ |Bl| := abs_nonneg Bl
  have hFB : -|Bl| ≤ D.faltingsHeight (D.cls E) := le_trans (neg_abs_le Bl) (hBl (D.cls E))
  obtain ⟨H, hHdef⟩ : ∃ H : ℝ, H = max 0 (23040 * 100 * (D.degOfDefinition E : ℝ)
      * (D.faltingsHeight (D.cls E) + C₈ * (D.degOfDefinition E : ℝ) ^ eps)) := ⟨_, rfl⟩
  have hHnn : (0:ℝ) ≤ H := by rw [hHdef]; exact le_max_left _ _
  have hXH : 23040 * 100 * (D.degOfDefinition E : ℝ)
      * (D.faltingsHeight (D.cls E) + C₈ * (D.degOfDefinition E : ℝ) ^ eps) ≤ H := by
    rw [hHdef]; exact le_max_right _ _
  have hFBd : (0:ℝ) ≤ (D.degOfDefinition E : ℝ) * D.faltingsHeight (D.cls E)
      + (D.degOfDefinition E : ℝ) * |Bl| := by nlinarith [hFB, hdpos]
  have hHle : H ≤ 23040 * 100 * ((D.degOfDefinition E : ℝ) * D.faltingsHeight (D.cls E))
      + 23040 * 100 * (C₈ * (D.degOfDefinition E : ℝ) ^ (1 + eps))
      + 23040 * 100 * ((D.degOfDefinition E : ℝ) * |Bl|) := by
    have hXeq : 23040 * 100 * (D.degOfDefinition E : ℝ)
          * (D.faltingsHeight (D.cls E) + C₈ * (D.degOfDefinition E : ℝ) ^ eps)
        = 23040 * 100 * ((D.degOfDefinition E : ℝ) * D.faltingsHeight (D.cls E))
          + 23040 * 100 * (C₈ * (D.degOfDefinition E : ℝ) ^ (1 + eps)) := by
      rw [← hdd]; ring
    have hdB : (0:ℝ) ≤ (D.degOfDefinition E : ℝ) * |Bl| := mul_nonneg hdpos.le hBnn
    have hCP : (0:ℝ) < C₈ * (D.degOfDefinition E : ℝ) ^ (1 + eps) :=
      mul_pos hC₈pos (Real.rpow_pos_of_pos hdpos _)
    rw [hHdef]
    rcases le_or_gt 0 (23040 * 100 * (D.degOfDefinition E : ℝ)
        * (D.faltingsHeight (D.cls E) + C₈ * (D.degOfDefinition E : ℝ) ^ eps)) with hX | hX
    · rw [max_eq_right hX, hXeq]; linarith
    · rw [max_eq_left hX.le]; linarith
  -- ★★`Lemma 4.2` を `h = 23040·d·deg∞` に適用する(原文 p.23)
  have hL42 := lemma_4_2 (D.multCard E) (D.multCard_pos E) (D.multPrime E) (D.multPrime_prime E)
      (D.localHt E) (D.localHt_pos E)
  have hsum := D.sum_localHt_eq E
  have hxbadle : (∑ p ∈ D.badPrimes E, Real.log p)
      ≤ 5/2 * (23040 * ((D.degOfDefinition E : ℝ) * D.degInf (D.cls E))) := by
    have h1 := hL42.1
    have h3 := hL42.2.2
    rw [hsum] at h1 h3
    have hb := D.sum_log_badPrimes_le E
    linarith
  have hAineq : D.degInf (D.cls E)
      ≤ 12 * (1 + 1/5) * D.faltingsHeight (D.cls E) + (A₁ + A₂) := by
    have h1 := hA₁ (D.cls E)
    have h2 := hA₂ (D.cls E)
    simp only at h1 h2
    linarith
  -- ★★★`Lemma 4.1` を `M = 1`･`1 + 6ϵ = 2` で適用する
  have hget : ∀ Bad : Finset ℕ, (∀ p ∈ Bad, Nat.Prime p) →
      ∃ l : ℕ, Nat.Prime l ∧ l ∉ S ∧ l ∉ Bad ∧ H ≤ (l:ℝ) ∧
        (l:ℝ) ≤ 2 * ((∑ p ∈ S, Real.log p) + (∑ p ∈ Bad, Real.log p)
          + (∑ p ∈ T, Real.log p)) + 8 * H := by
    intro Bad hBadp
    have hAllp : ∀ p ∈ S ∪ Bad ∪ T, Nat.Prime p := by
      intro p hp
      rcases Finset.mem_union.1 hp with hp' | hp'
      · rcases Finset.mem_union.1 hp' with hp'' | hp''
        · exact hSprime p hp''
        · exact hBadp p hp''
      · exact hTprime p hp'
    have hge : (∑ p ∈ T, Real.log p) ≤ ∑ p ∈ S ∪ Bad ∪ T, Real.log p :=
      Finset.sum_le_sum_of_subset_of_nonneg Finset.subset_union_right
        (fun p _ _ => Real.log_natCast_nonneg p)
    have hAx : xeps < ∑ p ∈ S ∪ Bad ∪ T, Real.log p := lt_of_lt_of_le hTgt hge
    obtain ⟨P, hcard, hPp, hPnot, hPb⟩ :=
      lemma_4_1 1 Nat.one_pos (1/6 : ℝ) xeps Ceps (by norm_num) hxepos hCepos (by norm_num) hCx
        (fun x hx _ => hi1 x hx) (fun x hx _ hxx => hi2 x hx hxx) hii H hHnn
        (S ∪ Bad ∪ T) hAllp hAx
    obtain ⟨a, ha⟩ := Finset.card_eq_one.1 hcard
    have hmem : a ∈ P := by rw [ha]; exact Finset.mem_singleton_self a
    have hnot := hPnot a hmem
    refine ⟨a, hPp a hmem, ?_, ?_, (hPb a hmem).1, ?_⟩
    · exact fun hS => hnot (Finset.mem_union_left _ (Finset.mem_union_left _ hS))
    · exact fun hBad => hnot (Finset.mem_union_left _ (Finset.mem_union_right _ hBad))
    · have h2 := (hPb a hmem).2
      have hnum : (1 + 6 * (1/6 : ℝ)) = 2 := by norm_num
      rw [hnum] at h2
      have hu1 := hunion (S ∪ Bad) T
      have hu2 := hunion S Bad
      linarith
  obtain ⟨lo, hlop, hloS, hlobad, hloH, hloub⟩ := hget (D.badPrimes E) (D.badPrimes_prime E)
  obtain ⟨lb, hlbp, hlbS, hlbram, hlbH, hlbub⟩ := hget (D.ramPrimes E) (D.ramPrimes_prime E)
  have hlbbad : lb ∉ D.badPrimes E := fun hc => hlbram (D.badPrimes_subset_ramPrimes E hc)
  obtain ⟨hlo1, hlo2⟩ := D.primeTo_badPrimes E lo hlop hlobad
  obtain ⟨hlb1, hlb2⟩ := D.primeTo_badPrimes E lb hlbp hlbbad
  have hlb3 : D.PrimeToRamification E lb := D.primeTo_ramPrimes E lb hlbp hlbram
  have hSL2lo : D.ImageContainsSL2 E lo :=
    h38 E lo hlop hExcE (Or.inl ⟨le_trans hXH hloH, hpot⟩)
  have hSL2lb : D.ImageContainsSL2 E lb :=
    h38 E lb hlbp hExcE (Or.inl ⟨le_trans hXH hlbH, hpot⟩)
  have hsurj : D.ImageSurjective E lb := D.imageSurjective_of_containsSL2 E lb hlbp hlb3 hSL2lb
  have hramle := D.sum_log_ramPrimes_le E
  have hclo := ABC3.Found.GenEll.cor4_numeric (D.degOfDefinition E : ℝ)
      (D.faltingsHeight (D.cls E)) (D.degInf (D.cls E)) (∑ p ∈ S, Real.log p)
      (∑ p ∈ D.badPrimes E, Real.log p) (∑ p ∈ T, Real.log p) 0 H (lo : ℝ)
      (A₁ + A₂) |Bl| C₈ ((D.degOfDefinition E : ℝ) ^ (1 + eps))
      hd1 hP1 hdP hAineq hxbadle hFB hBnn hC₈pos hxTnn hHle (by linarith [hloub])
  have hclb := ABC3.Found.GenEll.cor4_numeric (D.degOfDefinition E : ℝ)
      (D.faltingsHeight (D.cls E)) (D.degInf (D.cls E)) (∑ p ∈ S, Real.log p)
      (∑ p ∈ D.badPrimes E, Real.log p) (∑ p ∈ T, Real.log p)
      (3 * (D.degOfDefinition E : ℝ) * D.logDiffMell (D.cls E)) H (lb : ℝ)
      (A₁ + A₂) |Bl| C₈ ((D.degOfDefinition E : ℝ) ^ (1 + eps))
      hd1 hP1 hdP hAineq hxbadle hFB hBnn hC₈pos hxTnn hHle (by linarith [hlbub, hramle])
  exact ⟨lo, lb, hlop, hlbp, hloS, hlbS, ⟨hlo1, hlo2, hlb1, hlb2, hlb3⟩,
    ⟨hSL2lo, hsurj⟩, by linarith [hclo], by linarith [hclb]⟩

/-! ## Corollary 4.4 —— compactly bounded 部分集合の完全 Galois 作用(★北極星) -/

/-- **[GenEll] Corollary 4.4**(Full Galois Actions for Compactly Bounded Subsets)。

原文 (GenEll p.23):
> Corollary 4.4. (Full Galois Actions for Compactly Bounded Subsets)

`K_V` を compactly bounded subset とすると**定数 `C > 0` と Galois-finite な `Exc` が存在して**、
`[E_L] ∈ K_V` かつ `[E_L] ∉ Exc`、`L` が最小定義体なら、`S` に属さない素数 `l_∘, l_•` が存在して
(a)(b) は `Corollary 4.3` と同一、
(c) `l_∘ ≤ 23040·100d·ht^Falt + 2x_S + C·d`、
`l_• ≤ 23040·100d·ht^Falt + 6d·log-diff + 2x_S + C·d`。

★★**`Corollary 4.3` との差は (c) の 2 箇所だけである**(目視で 1 文字ずつ照合):
係数が `900d` → **`100d`**、末項が `C·d^{1+ϵ}` → **`C·d`**。
★すなわち **compactly bounded という仮定が `d` の指数を `1+ϵ` から `1` へ下げる**。
`ϵ` が statement から消えるのもこのためである。

★原文の証明は 3 行——「`Corollary 4.3` とまったく同様(むしろ易しい)。ただし
`Theorem 3.8` の条件 (a) の代わりに**条件 (b)** を使い、`Lemma 4.1` を使うときは
**`h` を `0` に取る**」。★`h = 0` にできることが `8h` の項を消して `C·d` を出す機構である。

★★**これが本トラックの北極星である**(`ResearchPaper/genell-goal.md`)。 -/
theorem cor_4_4 (D : EllModuliData) (KV : Set D.EllClass) (hKV : D.CompactlyBounded KV) :
    ∃ C : ℝ, 0 < C ∧ ∃ Exc : Set D.EllClass, D.GaloisFinite Exc ∧
      ∀ (E : D.Curve) (S : Finset ℕ), (∀ p ∈ S, p.Prime) →
        D.MinimalField E → D.cls E ∈ KV → D.cls E ∉ Exc →
        ∃ lo lb : ℕ, Nat.Prime lo ∧ Nat.Prime lb ∧ lo ∉ S ∧ lb ∉ S
          ∧ (D.PrimeToMultPrimes E lo ∧ D.PrimeToLocalHeights E lo
              ∧ D.PrimeToMultPrimes E lb ∧ D.PrimeToLocalHeights E lb
              ∧ D.PrimeToRamification E lb)
          ∧ (D.ImageContainsSL2 E lo ∧ D.ImageSurjective E lb)
          ∧ (lo : ℝ) ≤ 23040 * 100 * (D.degOfDefinition E : ℝ) * D.faltingsHeight (D.cls E)
              + 2 * (∑ p ∈ S, Real.log p) + C * (D.degOfDefinition E : ℝ)
          ∧ (lb : ℝ) ≤ 23040 * 100 * (D.degOfDefinition E : ℝ) * D.faltingsHeight (D.cls E)
              + 6 * (D.degOfDefinition E : ℝ) * D.logDiffMell (D.cls E)
              + 2 * (∑ p ∈ S, Real.log p) + C * (D.degOfDefinition E : ℝ) := by
  classical
  -- ★原文 p.23:『`Theorem 3.8` の条件 (a) の代わりに**条件 (b)** を使う』
  -- ——条件 (b) は `ϵ` を使わないので `ϵ := 1` でよい
  obtain ⟨C₈, hC₈pos, Exc, hExc, h38⟩ := theorem_3_8 D KV hKV 1 one_pos
  obtain ⟨hdegbd, hhtbd, -, -⟩ := prop_3_4 D (1/5 : ℝ) (by norm_num)
  obtain ⟨A₁, hA₁⟩ := hdegbd
  obtain ⟨A₂, hA₂⟩ := hhtbd
  obtain ⟨Bl, hBl⟩ := D.faltingsHeight_bddBelow
  obtain ⟨xeps, Ceps, hxepos, hCepos, hCx, hi1, hi2, hii⟩ :=
    ABC3.Found.GenEll.exists_cond_i_ii 1
  obtain ⟨T₀, hT₀prime, hT₀gt⟩ := ABC3.Found.GenEll.exists_finset_primes_sum_log_gt xeps
  -- ★★条件 (b) の『`30` と素』のために `2, 3, 5` を除外集合に入れておく
  obtain ⟨T, hTdef⟩ : ∃ T : Finset ℕ, T = T₀ ∪ ({2, 3, 5} : Finset ℕ) := ⟨_, rfl⟩
  have hTprime : ∀ p ∈ T, Nat.Prime p := by
    intro p hp
    rw [hTdef] at hp
    rcases Finset.mem_union.1 hp with hp' | hp'
    · exact hT₀prime p hp'
    · fin_cases hp' <;> norm_num
  have hT₀sub : T₀ ⊆ T := by rw [hTdef]; exact Finset.subset_union_left
  have hTgt : xeps < ∑ p ∈ T, Real.log p :=
    lt_of_lt_of_le hT₀gt (Finset.sum_le_sum_of_subset_of_nonneg hT₀sub
      (fun p _ _ => Real.log_natCast_nonneg p))
  have hxTnn : (0:ℝ) ≤ ∑ p ∈ T, Real.log p :=
    Finset.sum_nonneg fun p _ => Real.log_natCast_nonneg p
  have hm2 : (2:ℕ) ∈ T := by rw [hTdef]; exact Finset.mem_union_right _ (by decide)
  have hm3 : (3:ℕ) ∈ T := by rw [hTdef]; exact Finset.mem_union_right _ (by decide)
  have hm5 : (5:ℕ) ∈ T := by rw [hTdef]; exact Finset.mem_union_right _ (by decide)
  have hunion : ∀ U V : Finset ℕ,
      (∑ p ∈ U ∪ V, Real.log p) ≤ (∑ p ∈ U, Real.log p) + ∑ p ∈ V, Real.log p := by
    intro U V
    have hui : (∑ p ∈ U ∪ V, Real.log p) + ∑ p ∈ U ∩ V, Real.log p
        = (∑ p ∈ U, Real.log p) + ∑ p ∈ V, Real.log p := Finset.sum_union_inter
    have hnn : (0:ℝ) ≤ ∑ p ∈ U ∩ V, Real.log p :=
      Finset.sum_nonneg fun p _ => Real.log_natCast_nonneg p
    linarith
  refine ⟨23040 * (5 * |A₁ + A₂| + 28 * |Bl|) + 2 * (∑ p ∈ T, Real.log p) + 1,
    by nlinarith [abs_nonneg (A₁ + A₂), abs_nonneg Bl, hxTnn], Exc, hExc, ?_⟩
  intro E S hSprime _hmin hKVE hExcE
  have hd1 : (1:ℝ) ≤ (D.degOfDefinition E : ℝ) := by exact_mod_cast D.degOfDefinition_pos E
  have hBnn : (0:ℝ) ≤ |Bl| := abs_nonneg Bl
  have hFB : -|Bl| ≤ D.faltingsHeight (D.cls E) := le_trans (neg_abs_le Bl) (hBl (D.cls E))
  have hL42 := lemma_4_2 (D.multCard E) (D.multCard_pos E) (D.multPrime E) (D.multPrime_prime E)
      (D.localHt E) (D.localHt_pos E)
  have hsum := D.sum_localHt_eq E
  have hxbadle : (∑ p ∈ D.badPrimes E, Real.log p)
      ≤ 5/2 * (23040 * ((D.degOfDefinition E : ℝ) * D.degInf (D.cls E))) := by
    have h1 := hL42.1
    have h3 := hL42.2.2
    rw [hsum] at h1 h3
    have hb := D.sum_log_badPrimes_le E
    linarith
  have hAineq : D.degInf (D.cls E)
      ≤ 12 * (1 + 1/5) * D.faltingsHeight (D.cls E) + (A₁ + A₂) := by
    have h1 := hA₁ (D.cls E)
    have h2 := hA₂ (D.cls E)
    simp only at h1 h2
    linarith
  -- ★★★`Lemma 4.1` を `M = 1`･`h = 0` で適用する(原文 p.23)
  have hget : ∀ Bad : Finset ℕ, (∀ p ∈ Bad, Nat.Prime p) →
      ∃ l : ℕ, Nat.Prime l ∧ l ∉ S ∧ l ∉ Bad ∧ l ∉ T ∧
        (l:ℝ) ≤ 2 * ((∑ p ∈ S, Real.log p) + (∑ p ∈ Bad, Real.log p)
          + (∑ p ∈ T, Real.log p)) := by
    intro Bad hBadp
    have hAllp : ∀ p ∈ S ∪ Bad ∪ T, Nat.Prime p := by
      intro p hp
      rcases Finset.mem_union.1 hp with hp' | hp'
      · rcases Finset.mem_union.1 hp' with hp'' | hp''
        · exact hSprime p hp''
        · exact hBadp p hp''
      · exact hTprime p hp'
    have hge : (∑ p ∈ T, Real.log p) ≤ ∑ p ∈ S ∪ Bad ∪ T, Real.log p :=
      Finset.sum_le_sum_of_subset_of_nonneg Finset.subset_union_right
        (fun p _ _ => Real.log_natCast_nonneg p)
    have hAx : xeps < ∑ p ∈ S ∪ Bad ∪ T, Real.log p := lt_of_lt_of_le hTgt hge
    obtain ⟨P, hcard, hPp, hPnot, hPb⟩ :=
      lemma_4_1 1 Nat.one_pos (1/6 : ℝ) xeps Ceps (by norm_num) hxepos hCepos (by norm_num) hCx
        (fun x hx _ => hi1 x hx) (fun x hx _ hxx => hi2 x hx hxx) hii 0 le_rfl
        (S ∪ Bad ∪ T) hAllp hAx
    obtain ⟨a, ha⟩ := Finset.card_eq_one.1 hcard
    have hmem : a ∈ P := by rw [ha]; exact Finset.mem_singleton_self a
    have hnot := hPnot a hmem
    refine ⟨a, hPp a hmem, ?_, ?_, ?_, ?_⟩
    · exact fun hS => hnot (Finset.mem_union_left _ (Finset.mem_union_left _ hS))
    · exact fun hBad => hnot (Finset.mem_union_left _ (Finset.mem_union_right _ hBad))
    · exact fun hT => hnot (Finset.mem_union_right _ hT)
    · have h2 := (hPb a hmem).2
      have hnum : (1 + 6 * (1/6 : ℝ)) = 2 := by norm_num
      rw [hnum] at h2
      have hu1 := hunion (S ∪ Bad) T
      have hu2 := hunion S Bad
      linarith
  obtain ⟨lo, hlop, hloS, hlobad, hloT, hloub⟩ := hget (D.badPrimes E) (D.badPrimes_prime E)
  obtain ⟨lb, hlbp, hlbS, hlbram, hlbT, hlbub⟩ := hget (D.ramPrimes E) (D.ramPrimes_prime E)
  have hlbbad : lb ∉ D.badPrimes E := fun hc => hlbram (D.badPrimes_subset_ramPrimes E hc)
  obtain ⟨hlo1, hlo2⟩ := D.primeTo_badPrimes E lo hlop hlobad
  obtain ⟨hlb1, hlb2⟩ := D.primeTo_badPrimes E lb hlbp hlbbad
  have hlb3 : D.PrimeToRamification E lb := D.primeTo_ramPrimes E lb hlbp hlbram
  have hcoplo : Nat.Coprime lo 30 :=
    ABC3.Found.GenEll.coprime_thirty_of_prime hlop
      (fun h => hloT (h ▸ hm2)) (fun h => hloT (h ▸ hm3)) (fun h => hloT (h ▸ hm5))
  have hcoplb : Nat.Coprime lb 30 :=
    ABC3.Found.GenEll.coprime_thirty_of_prime hlbp
      (fun h => hlbT (h ▸ hm2)) (fun h => hlbT (h ▸ hm3)) (fun h => hlbT (h ▸ hm5))
  have hSL2lo : D.ImageContainsSL2 E lo := h38 E lo hlop hExcE (Or.inr ⟨hKVE, hlo2, hcoplo⟩)
  have hSL2lb : D.ImageContainsSL2 E lb := h38 E lb hlbp hExcE (Or.inr ⟨hKVE, hlb2, hcoplb⟩)
  have hsurj : D.ImageSurjective E lb := D.imageSurjective_of_containsSL2 E lb hlbp hlb3 hSL2lb
  have hramle := D.sum_log_ramPrimes_le E
  have hclo := ABC3.Found.GenEll.cor44_numeric (D.degOfDefinition E : ℝ)
      (D.faltingsHeight (D.cls E)) (D.degInf (D.cls E)) (∑ p ∈ S, Real.log p)
      (∑ p ∈ D.badPrimes E, Real.log p) (∑ p ∈ T, Real.log p) 0 (lo : ℝ)
      (A₁ + A₂) |Bl| (D.degOfDefinition E : ℝ)
      hd1 hd1 le_rfl hAineq hxbadle hFB hBnn hxTnn (by linarith [hloub])
  have hclb := ABC3.Found.GenEll.cor44_numeric (D.degOfDefinition E : ℝ)
      (D.faltingsHeight (D.cls E)) (D.degInf (D.cls E)) (∑ p ∈ S, Real.log p)
      (∑ p ∈ D.badPrimes E, Real.log p) (∑ p ∈ T, Real.log p)
      (3 * (D.degOfDefinition E : ℝ) * D.logDiffMell (D.cls E)) (lb : ℝ)
      (A₁ + A₂) |Bl| (D.degOfDefinition E : ℝ)
      hd1 hd1 le_rfl hAineq hxbadle hFB hBnn hxTnn (by linarith [hlbub, hramle])
  exact ⟨lo, lb, hlop, hlbp, hloS, hlbS, ⟨hlo1, hlo2, hlb1, hlb2, hlb3⟩,
    ⟨hSL2lo, hsurj⟩, by linarith [hclo], by linarith [hclb]⟩

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def lemma_4_1.src : Source :=
  { paper := "GenEll", pdfPage := 20, item := "Lemma 4.1",
    sectionId := "genell-lemma-4-1" }

/-- ★**空リストは省略ではなく主張である。**

原文 p.21 の証明は**背理法と 6 行の不等式計算**に尽きる。外部文献も別論文の項目も引かない
——条件 (i)(ii) が**仮説として与えられている**からである。
★最後の矛盾は `ϵ·x_A > ϵ·x_ϵ > C_ϵ` と `5(1+δ)/4 < 4` から出る(目視確認 2026-08-16)。

★ただし原文は末尾に「`y_A, (1+δ)h ∈ ℝ′_{>0}` を仮定してよい——`δ, h` をわずかに大きい
実数で置き換えればよく、議論に実質的な影響はない」と書いている。
★**この『実質的な影響はない』は段の飛躍であるが、`.needs` に書くべき外部依存ではない**
——同じ補題の証明の中の話なので、実装時に自分で埋める。 -/
def lemma_4_1.needs : List ProofObligation := []

def remark_4_1_1.src : Source :=
  { paper := "GenEll", pdfPage := 21, item := "Remark 4.1.1",
    sectionId := "genell-rem-4-1-1" }

def remark_4_1_1.needs : List ProofObligation :=
  [ .citation "[Edw]" "Riemann's Zeta Function, p. 76(素数定理の帰結として θ(x)/x → 1)"
      (.absent "Mathlib/NumberTheory/Chebyshev.lean は Chebyshev 型の評価(theta_ge' / psi_ge' / pi_le_log4_mul_div など)を持つが、θ(x)/x → 1 は無い(2026-08-16、同ファイルの全定理名を確認)。公開プロジェクト PrimeNumberTheoremAnd が持つが、使う前に clone して sorry を数えること") 21,
    .implicitStep
      "★第 1 項(条件 (ii) を満たす ϵ, x_ϵ, C_ϵ の存在)は『entirely elementary』とだけ書かれている。実際 M·log(x) ≤ ϵ·x は log の増大度から従い、mathlib の Real.isLittleO_log_id_atTop 相当で届く" 21 ]

def lemma_4_2.src : Source :=
  { paper := "GenEll", pdfPage := 21, item := "Lemma 4.2",
    sectionId := "genell-lemma-4-2" }

/-- ★**空リストは省略ではなく主張である。** 原文の証明は 2 行で、外部依存を持たない。 -/
def lemma_4_2.needs : List ProofObligation := []

def cor_4_3.src : Source :=
  { paper := "GenEll", pdfPage := 22, item := "Corollary 4.3",
    sectionId := "genell-cor-4-3" }

def cor_4_3.needs : List ProofObligation :=
  [ .otherPaper "[GenEll]" "Theorem 3.8(条件 (a) を使う)" 19,
    .otherPaper "[GenEll]" "Lemma 4.1(指定された大きさの素数の存在)" 20,
    .otherPaper "[GenEll]" "Lemma 4.2(初等的な評価)" 21,
    .otherPaper "[GenEll]" "Proposition 3.4(deg_∞ と ht^Falt の比較)" 17,
    .otherPaper "[GenEll]" "Definition 1.5, (iii)(log-diff_{M̄_ell})" 8,
    .folklore
      "原文「the well-known fact that the field extension ℚ(ζ_{l^∞})/ℚ obtained by adjoining the l-th power roots of unity to ℚ is totally ramified over the prime l, hence linearly disjoint from the extension L/ℚ」——これで『像が SL₂ を含む』⇔『全射』が出る" 22,
    .implicitStep
      "★原文は『the primes appearing in the arithmetic divisor that gives rise to log-diff_{M̄_ell} appear with multiplicity ≥ one less than the ramification indices of L/ℚ』を『as is easily verified, by considering the trace of an extension of number fields』で済ませている" 23 ]

def cor_4_4.src : Source :=
  { paper := "GenEll", pdfPage := 23, item := "Corollary 4.4",
    sectionId := "genell-cor-4-4" }

/-- ★原文の証明は 3 行で、`Corollary 4.3` との差だけを述べている。 -/
def cor_4_4.needs : List ProofObligation :=
  [ .otherPaper "[GenEll]" "Theorem 3.8(★条件 (b) を使う——Corollary 4.3 は条件 (a))" 19,
    .otherPaper "[GenEll]" "Corollary 4.3(証明は『entirely similar to [but slightly easier than]』)" 22,
    .otherPaper "[GenEll]" "Lemma 4.1(★h を 0 に取る——これが C·d^{1+ϵ} を C·d に落とす機構)" 20,
    .otherPaper "[GenEll]" "Example 1.3, (ii)(compactly bounded subset)" 5,
    .implicitStep
      "★原文は 3 行で『Corollary 4.3 とまったく同様』と述べるだけで、係数が 900d から 100d に、末項が C·d^{1+ϵ} から C·d に変わる理由を明示的には追っていない。**この差分の検証が実装時の実質である**" 23 ]

end ABC3.Skeleton.GenEll
