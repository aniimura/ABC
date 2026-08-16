import ABC3.Meta.Claim
import ABC3.Interface.GenEll.EllModuli
import ABC3.Found.GenEll.Elementary
import ABC3.Found.GenEll.PrimesOfSize
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
第 1 項(初等的な部分)は今すぐ証明できる。
★**第 2 項が素数定理そのもの**であり、mathlib に無い(上の docstring で実測)。 -/
theorem remark_4_1_1 :
    (∀ (M : ℕ) (eps : ℝ), 0 < eps →
        ∃ xeps : ℝ, 0 < xeps ∧ ∀ x : ℝ, xeps ≤ x → (M : ℝ) * Real.log x ≤ eps * x)
  ∧ Filter.Tendsto (fun x : ℝ => Chebyshev.theta x / x) Filter.atTop (nhds 1) := by
  sorry

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
  sorry

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
  sorry

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
