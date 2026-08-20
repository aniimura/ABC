/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop56Sec
import ABC3.Found.FrdI.Prop56MpAbs

/-!
# [FrdI] Proposition 5.6 —— `u` の構成に向けた部品

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.106–107。

原文 (FrdI p.107):
> product” of the up].

原文 (FrdI p.107):
> — which [by the total epimorphicity of C] implies that wl

## ★★原文の `u` の構成(4 段)

`prop_5_6_claim`(`Prop56Sec.lean`)により、`Proposition 5.6` は
**`u ∈ 𝒪^×(A)` で `u · φ_p · u⁻¹ = φ'_p`(すべての素数 `p`)なるものの存在**に
還元されている。原文はそれを次の 4 段で作る:

1. `φ'_l ≈_p v_l · φ_l`(`v_l ∈ 𝒪^×(A)`)—— **本ファイルの `exists_otimes_of_frobType_same_deg`**
2. `w · φ_l · w⁻¹ ≈_p w^{1-l} · φ_l` —— **本ファイルの `conj_frobNorm`**
3. `u_p` の存在(`u^{1-p} = v_p` を pro-`p` 群で解く)——
   `Found/ProL/DivByP.lean` の `exists_pow_one_sub`
4. `w_l ≈_p w_l^p ⟹ w_l ≈_p 1` —— 同 `eq_one_of_eq_pow_self`、
   および `∏_p 𝒪^×(A)[p] ≅ 𝒪^×(A)`(`decompEquiv`)で無限積を取る

★★本ファイルは **1 と 2** を置く。3 と 4 の群論側は済んでいるので、
残るのは **`E_{M_p}` の上での計算と無限積の組み立て**(鎖 `prop56` の `p56-limit`)である。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

/-- ★★★**同じ次数の base-identity な Frobenius 型自己射は、単元だけずれる**。

原文 (FrdI p.107):
> product” of the up].

★★これが原文の「`φ'_l ≈_p v_l · φ_l`」の中身である。
★`Definition 1.3, (ii)` の**本質的一意性**(`frobDegUniq`)から出る ——
そこで得た `β` は底が恒等で同型なので `𝒪^×(A)` に入る。 -/
theorem exists_otimes_of_frobType_same_deg (G : Frobenioid P) {A : C} (φ ψ : End A)
    (hφb : IsBaseIdentity P ((φ : A ⟶ A))) (hψb : IsBaseIdentity P ((ψ : A ⟶ A)))
    (hφ : IsFrobeniusType P ((φ : A ⟶ A))) (hψ : IsFrobeniusType P ((ψ : A ⟶ A)))
    (hd : P.degFr ((φ : A ⟶ A)) = P.degFr ((ψ : A ⟶ A))) :
    ∃ v, v ∈ OTimes P A ∧ ψ = v * φ := by
  obtain ⟨β, hβiso, hβ⟩ := G.core.frobDegUniq A A A (φ : A ⟶ A) (ψ : A ⟶ A) hφ hψ hd
  haveI := hβiso
  have hbβ : P.Base β = P.Base (𝟙 A) := by
    have h := congrArg P.Base hβ
    rw [P.Base_comp, hφb, hψb, P.Base_id, Category.id_comp] at h
    rw [h, P.Base_id]
  refine ⟨(β : End A), ⟨⟨hbβ, degFr_of_isIso P β⟩, (isUnit_iff_isIso (β : End A)).mpr hβiso⟩, ?_⟩
  exact hβ.symm

/-- ★★★**原文 p.107 の「`w · φ_l · w⁻¹ ≈ w^{1-l} · φ_l`」**。

原文 (FrdI p.107):
> Since, for w ∈O×(A)[p], we have, for l ∈Primes, w·φl ·w−1

★Frobenius-normalized(`φ · α = α^{deg φ} · φ`)を `α := w⁻¹` に当て、
結合則で `w · (φ · w⁻¹) = w · (w⁻¹)^l · φ` にするだけである。 -/
theorem conj_frobNorm {A : C} (hfn : IsFrobeniusNormalized P A) {φ w wi : End A}
    (hφ : IsBaseIdentity P ((φ : A ⟶ A))) (hwi : wi ∈ OTri P A) :
    w * φ * wi = w * wi ^ ((P.degFr ((φ : A ⟶ A)) : ℕ+) : ℕ) * φ := by
  rw [mul_assoc, frobNorm_end hfn hφ hwi, ← mul_assoc]

/-! ### ★出典の紐付け -/

/-- ★locator —— `Proposition 5.6` の `φ'_l = v_l · φ_l`。 -/
def exists_otimes_of_frobType_same_deg.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 107,
    item := "Proposition 5.6 — φ'_l は φ_l と単元だけずれる",
    sectionId := "frdi-prop-5-6" }

/-- ★locator —— `Proposition 5.6` の `w · φ_l · w⁻¹ = w^{1-l} · φ_l`。 -/
def conj_frobNorm.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 107,
    item := "Proposition 5.6 — w·φ_l·w⁻¹ = w^{1-l}·φ_l",
    sectionId := "frdi-prop-5-6" }

end ABC3.Found.FrdI
