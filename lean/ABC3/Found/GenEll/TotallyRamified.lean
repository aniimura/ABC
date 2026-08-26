import ABC3.Found.GenEll.DifferentKummer
import Mathlib.RingTheory.Nakayama
import Mathlib.FieldTheory.Minpoly.IsIntegrallyClosed

/-!
# [GenEll] Proposition 1.7 —— 全分岐の構造定理(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.10。

原文 (GenEll p.10):
> Σ” of log-condE, log-condD is ≈0 [cf. Remark 1.5.1], while [again by the elementary

## ★★★★★★★★これが段 1 の最後の穴だった

`Proposition 1.7, (i)` の証明が畳んでいる "elementary claim" の段 1 は
**全分岐拡大の構造定理**である:

    B が A の全分岐拡大で λ が B の素元なら、B = A[λ] であり
    minpoly A λ は次数 e の Eisenstein 多項式である

`DifferentKummer.lean` の `mem_differentIdeal_of_totallyRamified_tame` は
これを**仮定として**受け取る形で書いてあった(`hmono`・`hdeg`)。ここでそれを埋める。

## ★★★実測(2026-08-27)—— mathlib に「全分岐」は無い

★**`totally ramified` という語が mathlib に 1 件も無い**(全文検索 0 件)。
`IsTotallyRamified` も無い。分岐の仮説から `Algebra.adjoin A {λ} = ⊤` を出す定理も無い。

| 要るもの | mathlib |
|---|---|
| `e·f = [L:K]` の局所版 | ★ある(`Ideal.ramificationIdx_mul_inertiaDeg_of_isLocalRing`) |
| Nakayama | ★ある(`Submodule.le_of_le_smul_of_le_jacobson_bot`) |
| `B = A[λ]` | ❌ 無い —— **本ファイルで自作する** |
| `minpoly A λ` の次数が `e` | ❌ 無い —— **本ファイルで自作する** |

★★唯一の `Algebra.adjoin R {β} = ⊤` は `IsLocalRing.exists_adjoin_eq_top` で、
`FormallyUnramified` すなわち**逆の場合**である。
その Nakayama 論法は `map_maximalIdeal` を消費するが、全分岐では偽(そこでは `𝔪_B^e`)。

## ★★★★★仮定を「初等的な形」で置く

『剰余体が一致する』を **`∀ b : B, ∃ a : A, b − a ∈ 𝔪_B`** と書いた。
★これは `inertiaDeg = 1` の初等版であり、Dedekind 環の道具立てを持ち込まずに済む。
★★`e = [L:K]` から `f = 1` を出すのは
`Ideal.ramificationIdx_mul_inertiaDeg_of_isLocalRing` の仕事で、**別の段**である。
-/

namespace ABC3.Found.GenEll

open Polynomial

/-! ## ★★★★★★★★`B = A[λ]` —— Nakayama -/

/-- ★★★★★★★★**全分岐なら素元 1 つで生成される**。

原文 (GenEll p.10):
> Σ” of log-condE, log-condD is ≈0 [cf. Remark 1.5.1], while [again by the elementary

仮定は 3 つ:
* `hres` —— 剰余体が一致する(`∀ b, ∃ a ∈ A, b − a ∈ 𝔪_B`)
* `hmax` —— `𝔪_B = (λ)`(`λ` が素元)
* `hpow` —— `λ^e` が `π` の単元倍(全分岐)

★証明は **Nakayama 1 回**である。`hres` と `hmax` を `e` 回繰り返すと
`b ∈ A[λ] + λ^e·B` が出て、`hpow` が `λ^e·B ⊆ 𝔪_A·B` を与える。
★★`Module.Finite A B` が `⊤` の有限生成性を供給し、
`IsLocalRing A` が `𝔪_A ≤ jacobson ⊥` を供給する。 -/
theorem adjoin_eq_top_of_residue_surj
    {A B : Type*} [CommRing A] [IsLocalRing A] [CommRing B] [IsLocalRing B]
    [Algebra A B] [Module.Finite A B]
    (lam : B) (e : ℕ) (pi : A)
    (hres : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ IsLocalRing.maximalIdeal B)
    (hmax : IsLocalRing.maximalIdeal B = Ideal.span {lam})
    (hpiA : pi ∈ IsLocalRing.maximalIdeal A)
    (hpow : ∃ u : B, IsUnit u ∧ lam ^ e = algebraMap A B pi * u) :
    Algebra.adjoin A {lam} = ⊤ := by
  classical
  set N : Submodule A B := Subalgebra.toSubmodule (Algebra.adjoin A {lam}) with hN
  have hpowN : ∀ n : ℕ, lam ^ n ∈ N := by
    intro n
    rw [hN]
    exact Subalgebra.pow_mem _ (Algebra.subset_adjoin (Set.mem_singleton lam)) n
  -- ★1 段: `b = a + λ·b′`
  have hstep : ∀ b : B, ∃ (a : A) (b' : B), b = algebraMap A B a + lam * b' := by
    intro b
    obtain ⟨a, ha⟩ := hres b
    rw [hmax, Ideal.mem_span_singleton] at ha
    obtain ⟨c, hc⟩ := ha
    exact ⟨a, c, by linear_combination hc⟩
  -- ★`n` 段: `b ∈ A[λ] + λ^n·B`
  have hiter : ∀ n : ℕ, ∀ b : B, ∃ z ∈ N, ∃ b' : B, b = z + lam ^ n * b' := by
    intro n
    induction n with
    | zero => intro b; exact ⟨0, N.zero_mem, b, by simp⟩
    | succ n ih =>
      intro b
      obtain ⟨z, hzN, b', hb⟩ := ih b
      obtain ⟨a, b'', hb'⟩ := hstep b'
      refine ⟨z + a • lam ^ n, N.add_mem hzN (N.smul_mem a (hpowN n)), b'', ?_⟩
      rw [hb, hb', Algebra.smul_def]
      ring
  -- ★`e` 段目で `λ^e·B ⊆ 𝔪_A·B` になる
  have htop : (⊤ : Submodule A B) ≤ N ⊔ (IsLocalRing.maximalIdeal A) • (⊤ : Submodule A B) := by
    intro b _
    obtain ⟨z, hzN, b', hb⟩ := hiter e b
    obtain ⟨u, _, hpe⟩ := hpow
    have hrw : lam ^ e * b' = pi • (u * b') := by
      rw [hpe, Algebra.smul_def]; ring
    rw [hb, hrw]
    exact Submodule.add_mem_sup hzN (Submodule.smul_mem_smul hpiA Submodule.mem_top)
  have hfg : (⊤ : Submodule A B).FG := Module.Finite.fg_top
  have hjac : IsLocalRing.maximalIdeal A ≤ (⊥ : Ideal A).jacobson := by
    rw [IsLocalRing.jacobson_eq_maximalIdeal ⊥ bot_ne_top]
  have hle := Submodule.le_of_le_smul_of_le_jacobson_bot hfg hjac htop
  refine Algebra.eq_top_iff.2 (fun x => ?_)
  have hx := hle (Submodule.mem_top (x := x))
  rwa [hN] at hx

/-! ## ★★★★★★`minpoly` の次数は拡大次数 -/

/-- ★★★★**`L = K(x)` なら `deg (minpoly K x) = [L:K]`**。 -/
theorem natDegree_minpoly_eq_finrank {K L : Type*} [Field K] [Field L] [Algebra K L]
    [FiniteDimensional K L] (x : L) (hgen : Algebra.adjoin K {x} = ⊤) :
    (minpoly K x).natDegree = Module.finrank K L := by
  have hint : IsIntegral K x := IsIntegral.of_finite K x
  have htop : (IntermediateField.adjoin K {x} : IntermediateField K L) = ⊤ := by
    refine IntermediateField.toSubalgebra_injective ?_
    rw [IntermediateField.adjoin_simple_toSubalgebra_of_integral hint, hgen]
    rfl
  rw [← IntermediateField.adjoin.finrank hint, htop]
  exact (IntermediateField.finrank_top' (F := K) (E := L)).symm ▸ rfl

/-- ★★★★★★**`deg (minpoly A λ) = [L:K]`** —— 整閉なら `A` 側と `K` 側で次数が変わらない。

★`minpoly.isIntegrallyClosed_eq_field_fractions` で `minpoly K (λ) = (minpoly A λ).map`
になり、`algebraMap A K` は単射なので次数は保たれる。

★★これが `mem_differentIdeal_of_totallyRamified_tame` の `hdeg` を供給する。 -/
theorem natDegree_minpoly_eq_finrank_of_adjoin
    (A K L : Type*) {B : Type*} [CommRing A] [Field K] [CommRing B] [Field L]
    [IsDomain A] [IsDomain B]
    [Algebra A K] [Algebra B L] [Algebra A B] [Algebra K L] [Algebra A L]
    [IsScalarTower A K L] [IsScalarTower A B L] [IsFractionRing A K]
    [FiniteDimensional K L] [IsIntegrallyClosed A]
    (lam : B) (hint : IsIntegral A lam)
    (hgen : Algebra.adjoin K {(algebraMap B L) lam} = ⊤) :
    (minpoly A lam).natDegree = Module.finrank K L := by
  have h1 := natDegree_minpoly_eq_finrank ((algebraMap B L) lam) hgen
  rw [minpoly.isIntegrallyClosed_eq_field_fractions K L hint] at h1
  rw [← h1]
  exact (Polynomial.natDegree_map_eq_of_injective (IsFractionRing.injective A K) _).symm

/-! ## ★★★★★★★★★★段 1 が閉じた —— 全分岐の仮説から `p ∈ 𝔡` へ -/

/-- ★★★★★★★★★★**[GenEll] Proposition 1.7, (i) の elementary claim、段 1(馴の側)**。

原文 (GenEll p.10):
> Σ” of log-condE, log-condD is ≈0 [cf. Remark 1.5.1], while [again by the elementary

`DifferentKummer.lean` の `mem_differentIdeal_of_totallyRamified_tame` は
`hmono`(`B = A[λ]`)と `hdeg`(`deg minpoly = e`)を**仮定として**受け取っていた。
★本定理はその 2 つを**全分岐の初等的な仮説から導いて**渡す。

仮説は 5 つで、いずれも「全分岐かつ馴」の初等的な言い換えである:

| 仮説 | 意味 |
|---|---|
| `hres` | 剰余体が一致する(`f = 1`) |
| `hmax` | `𝔪_B = (λ)`(`λ` が素元) |
| `hpow` | `λ^e` が `π` の単元倍(`e = ramificationIdx`) |
| `hfinrank` | `[L:K] = e`(全分岐) |
| `hunit` | `e` が `B` で単元(馴) |

★★これで `Proposition 1.7` の elementary claim の 6 段のうち**段 1 が閉じた**。 -/
theorem mem_differentIdeal_of_totallyRamified_tame'
    (A : Type*) (K : Type*) (L : Type*) {B : Type*} [CommRing A] [Field K] [CommRing B] [Field L]
    [Algebra A K] [Algebra B L] [Algebra A B] [Algebra K L] [Algebra A L]
    [IsScalarTower A K L] [IsScalarTower A B L] [IsDomain A] [IsFractionRing A K]
    [FiniteDimensional K L] [Algebra.IsSeparable K L] [IsIntegralClosure B A L]
    [IsIntegrallyClosed A] [IsDedekindDomain B] [Module.IsTorsionFree A B]
    [IsLocalRing A] [IsLocalRing B] [Module.Finite A B]
    (e m : ℕ) (he : 0 < e) (lam pp : B) (pi : A) (hlam0 : lam ≠ 0)
    (hres : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ IsLocalRing.maximalIdeal B)
    (hmax : IsLocalRing.maximalIdeal B = Ideal.span {lam})
    (hpiA : pi ∈ IsLocalRing.maximalIdeal A)
    (hpow : ∃ u : B, IsUnit u ∧ lam ^ e = algebraMap A B pi * u)
    (hdvdA : ∀ x : A, ¬ IsUnit x → lam ^ e ∣ algebraMap A B x)
    (hloc : ∀ x : A, ¬ IsUnit (algebraMap A B x) → ¬ IsUnit x)
    (hint : IsIntegral A lam)
    (hfinrank : Module.finrank K L = e)
    (hunit : IsUnit (e : B))
    (hgen : Algebra.adjoin K {(algebraMap B L) lam} = ⊤)
    (hle : e - 1 ≤ m) (hdvdp : lam ^ m ∣ pp) :
    pp ∈ differentIdeal A B := by
  have hlamm : lam ∈ IsLocalRing.maximalIdeal B := by
    rw [hmax]
    exact Ideal.mem_span_singleton_self lam
  have hmono : Algebra.adjoin A {lam} = ⊤ :=
    adjoin_eq_top_of_residue_surj lam e pi hres hmax hpiA hpow
  have hdeg : (minpoly A lam).natDegree = e := by
    rw [natDegree_minpoly_eq_finrank_of_adjoin A K L lam hint hgen, hfinrank]
  exact mem_differentIdeal_of_totallyRamified_tame A K L e m he lam pp hlam0 hlamm
    hdvdA hloc hint hdeg hunit hmono hgen hle hdvdp

/-! ## ★出典の紐付け(`.src`) -/

def adjoin_eq_top_of_residue_surj.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7, (i) の elementary claim(全分岐なら B = A[λ])",
    sectionId := "genell-prop-1-7" }

def natDegree_minpoly_eq_finrank_of_adjoin.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7, (i) の elementary claim(minpoly A λ の次数が e)",
    sectionId := "genell-prop-1-7" }

def mem_differentIdeal_of_totallyRamified_tame'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7, (i) の elementary claim(段 1 の馴の側を全分岐の初等的な仮説から)",
    sectionId := "genell-prop-1-7" }

end ABC3.Found.GenEll
