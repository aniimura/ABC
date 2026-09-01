/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import Mathlib.FieldTheory.PrimitiveElement
import Mathlib.FieldTheory.Normal.Basic
import Mathlib.FieldTheory.IntermediateField.Adjoin.Algebra
import Mathlib.FieldTheory.IsAlgClosed.AlgebraicClosure
import ABC3.Meta.Claim

/-!
# 第 1214 ブロック —— **共役の個数で次数を上から抑える**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か——`[M:L] < l` へ向かう段

第 1209 で「残る義務は `[M:L] < l` の 1 つ」になった。
☆原文の道では `M = L(H)` の Galois 群が `𝔽_l^×` に埋まるので
`[M:L] ∣ l−1 < l` である。

★その議論の**数え上げの部分**を、`Galois` を使わない形で取る:

    `M →ₐ[L] L̄` が有限集合 `s` に単射で入る  ⟹  `[M:L] ≤ #s`

☆標数 0（分離的）なら `#(M →ₐ[L] L̄) = [M:L]`（mathlib の `AlgHom.card`）。

★★★あとは「`φ : M →ₐ[L] L̄` は `Q ↦ c·Q` の `c ∈ 𝔽_l^×` で決まる」を示せば
`[M:L] ≤ l−1 < l` が出る——第 1205 と第 1213 がその材料である。
-/

namespace ABC3.Found.GenEll

open ABC3.Meta

/-- ★★★★★★★★★★★★
**共役の個数で次数を上から抑える**——★**無条件**（第 1214）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆標数 0（分離的）なら `#(M →ₐ[L] L̄) = [M:L]` なので、
`M →ₐ[L] L̄` が有限集合 `s` に単射で入れば `[M:L] ≤ #s` である。

★★★これが `[M:L] < l`（第 1209 が残した唯一の義務）へ向かう
**数え上げの段**である。 -/
theorem finrank_le_of_algHom_mem_finset
    {L M Lbar : Type*} [Field L] [Field M] [Field Lbar]
    [Algebra L M] [Algebra L Lbar] [IsAlgClosed Lbar]
    [FiniteDimensional L M] [Algebra.IsSeparable L M]
    {α : Type*} [DecidableEq α] (f : (M →ₐ[L] Lbar) → α)
    (hinj : Function.Injective f) (s : Finset α) (hs : ∀ φ, f φ ∈ s) :
    Module.finrank L M ≤ s.card := by
  classical
  have hcard : Fintype.card (M →ₐ[L] Lbar) = Module.finrank L M := AlgHom.card L M Lbar
  rw [← hcard, ← Finset.card_univ]
  exact Finset.card_le_card_of_injOn f (fun a _ => hs a) (fun a _ b _ h => hinj h)

/-! ## ★★★★★★★★★★★★中間体の埋め込みは自己同型に延びる -/

/-- ★★★★★★★★★★★★
**中間体からの埋め込みは大域の自己同型に延びる**——★**無条件**（第 1215）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆mathlib の `AlgHom.liftNormal`（`E/Kᵢ/F` の塔で `E/F` が正規なら
`ϕ : K₁ →ₐ[F] K₂` は `E →ₐ[F] E` に延びる）と
`AlgHom.normal_bijective`（正規拡大の自己準同型は全単射）を繋いだもの。

★★★これで「`φ : M →ₐ[L] L̄` は `Gal(L̄/L)` の元から来る」が言えるので、
第 1205（`σ Q = c • Q`）が `φ` にも当たる
——第 1214 の数え上げと合わせて `[M:L] ≤ l−1` へ向かう。 -/
theorem exists_algEquiv_extend {L Lbar M : Type*} [Field L] [Field Lbar] [Field M]
    [Algebra L Lbar] [Normal L Lbar]
    [Algebra L M] [Algebra M Lbar] [IsScalarTower L M Lbar]
    (φ : M →ₐ[L] Lbar) :
    ∃ σ : Lbar ≃ₐ[L] Lbar, ∀ x : M, σ (algebraMap M Lbar x) = φ x := by
  refine ⟨AlgEquiv.ofBijective (φ.liftNormal Lbar)
    (AlgHom.normal_bijective L Lbar Lbar _), fun x => ?_⟩
  show φ.liftNormal Lbar (algebraMap M Lbar x) = φ x
  rw [AlgHom.liftNormal_commutes]
  rfl

/-! ## ★★★★★★★★★★★★★★★★不変量が生成元の上で決めるなら次数が抑えられる -/

/-- ★★★★★★★★★★★★★★★★
**不変量が生成元の上で `φ` を決めるなら次数が抑えられる**——★**無条件**（第 1216）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`g` が `φ` を生成集合 `s` の上で決めるなら、`AlgHom.ext_of_adjoin_eq_top` で
`g` は単射であり、第 1214 で `[M:L] ≤ #t` が出る。

★★★これが `[M:L] ≤ l−1` の**最終形**である
——`g φ ≔ c(φ)`（`σ Q = c • Q` の `c`、第 1205・第 1215）、
`t ≔ mod l の 0 でない元` と取ればよい。 -/
theorem finrank_le_of_determines_on_generators
    {L M Lbar : Type*} [Field L] [Field M] [Field Lbar]
    [Algebra L M] [Algebra L Lbar] [IsAlgClosed Lbar]
    [FiniteDimensional L M] [Algebra.IsSeparable L M]
    {α : Type*} [DecidableEq α] (s : Set M) (hs : Algebra.adjoin L s = ⊤)
    (g : (M →ₐ[L] Lbar) → α)
    (hg : ∀ φ₁ φ₂ : M →ₐ[L] Lbar, g φ₁ = g φ₂ → Set.EqOn φ₁ φ₂ s)
    (t : Finset α) (ht : ∀ φ, g φ ∈ t) :
    Module.finrank L M ≤ t.card :=
  finrank_le_of_algHom_mem_finset g
    (fun φ₁ φ₂ h => AlgHom.ext_of_adjoin_eq_top hs (hg φ₁ φ₂ h)) t ht

/-! ## ★★★★★★★★★★★★生成集合を `Algebra.adjoin` の形で出す -/

/-- ★★★★★★★★★★★★
**`L(T)` の中で `T` は `Algebra.adjoin` としても生成する**——★**無条件**（第 1217）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`T` が代数的なら `(IntermediateField.adjoin L T).toSubalgebra = Algebra.adjoin L T`
（mathlib の `adjoin_toSubalgebra_of_isAlgebraic`）なので、
`M ≔ L(T)` の中で `T` の引き戻しは `Algebra.adjoin` で `⊤` を生む。

★★★これが第 1216 が受け取る `hs : Algebra.adjoin L s = ⊤` を**実際に取る**段である。 -/
theorem algebra_adjoin_preimage_eq_top {L Lbar : Type*} [Field L] [Field Lbar]
    [Algebra L Lbar] (T : Set Lbar) (hT : ∀ x ∈ T, IsAlgebraic L x) :
    Algebra.adjoin L
        (((↑) : ↥(IntermediateField.adjoin L T) → Lbar) ⁻¹' T) = ⊤ := by
  set M := IntermediateField.adjoin L T with hMdef
  set s : Set M := ((↑) : M → Lbar) ⁻¹' T with hsdef
  have himg : (M.val : M →ₐ[L] Lbar) '' s = T := by
    ext y
    constructor
    · rintro ⟨z, hz, rfl⟩
      exact hz
    · intro hy
      exact ⟨⟨y, IntermediateField.subset_adjoin L T hy⟩, hy, rfl⟩
  have hmap : Subalgebra.map (M.val : M →ₐ[L] Lbar) (Algebra.adjoin L s)
      = Algebra.adjoin L T := by
    rw [AlgHom.map_adjoin, himg]
  have htop : Subalgebra.map (M.val : M →ₐ[L] Lbar) ⊤ = Algebra.adjoin L T := by
    rw [Algebra.map_top, IntermediateField.range_val,
      IntermediateField.adjoin_toSubalgebra_of_isAlgebraic hT]
  exact Subalgebra.map_injective (M.val).injective (hmap.trans htop.symm)

/-! ## ★出典の紐付け(`.src`) -/

def algebra_adjoin_preimage_eq_top.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(L(T) の中で T は Algebra.adjoin としても生成する。★無条件)",
    sectionId := "genell-lemma-3-5" }

def finrank_le_of_determines_on_generators.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(不変量が生成元の上で φ を決めるなら次数が抑えられる。★無条件)",
    sectionId := "genell-lemma-3-5" }

def finrank_le_of_determines_on_generators.needs : List ProofObligation :=
  [ .implicitStep
      ("★★★★**2026-09-02（第 1216）**——`[M:L] ≤ l−1` の**最終形**である。" ++
       "☆`g φ ≔ c(φ)`（`σ Q = c • Q` の `c`、第 1205・第 1215）、" ++
       "`t ≔ mod l の 0 でない元` と取ればよい。" ++
       "★残るのは `M` の生成集合を `Algebra.adjoin` の形で用意する段である。") 2 ]

def exists_algEquiv_extend.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(中間体からの埋め込みは大域の自己同型に延びる。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_algEquiv_extend.needs : List ProofObligation :=
  [ .implicitStep
      ("★★★★**2026-09-02（第 1215）**——「`φ : M →ₐ[L] L̄` は `Gal(L̄/L)` の" ++
       "元から来る」の段である。☆mathlib の `AlgHom.liftNormal` と" ++
       "`AlgHom.normal_bijective` を繋いだ。" ++
       "★これで第 1205（`σ Q = c • Q`）が `φ` にも当たるので、" ++
       "第 1214 の数え上げと合わせて `[M:L] ≤ l−1` へ向かう。") 2 ]

def finrank_le_of_algHom_mem_finset.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(共役の個数で次数を上から抑える。★無条件)",
    sectionId := "genell-lemma-3-5" }

def finrank_le_of_algHom_mem_finset.needs : List ProofObligation :=
  [ .implicitStep
      ("★★★★**2026-09-02（第 1214）**——第 1209 が残した唯一の義務 `[M:L] < l` へ" ++
       "向かう**数え上げの段**である。☆標数 0（分離的）なら" ++
       "`#(M →ₐ[L] L̄) = [M:L]`（mathlib の `AlgHom.card`）。" ++
       "★あとは「`φ : M →ₐ[L] L̄` は `Q ↦ c·Q` の `c ∈ 𝔽_l^×` で決まる」を" ++
       "示せば `[M:L] ≤ l−1 < l` が出る——第 1205 と第 1213 がその材料である。") 3 ]

end ABC3.Found.GenEll
