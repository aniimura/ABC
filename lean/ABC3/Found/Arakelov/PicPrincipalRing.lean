import ABC3.Found.Arakelov.PicIdealMul

/-!
# Arakelov (B2) 第 191 ブロック —— **可逆な単項イデアルは自由**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★主因子の 4 欄に要る環レベルの補題

    Module.Invertible R I  ∧  I 単項  ⟺  I ≃ₗ[R] R

★★「単項 ⟹ 自由」に**可逆性が要る**のが要点である。`R = k[x]/(x²)`、`I = (x)` は
単項だが `I ≅ k`(1 次元)で `R`(2 次元)とは同型でない。★可逆性(= 生成元が
非零因子)を課すと消える。**Interface の 3 件目の誤り**(§9-209)はここであった。

## ★★mathlib の `bijective_of_surjective` が効いた

    Module.Invertible.bijective_of_surjective :
      [Module.Invertible R N] → Surjective (f : M →ₗ[R] N) → Bijective f

★`R → I`(`r ↦ r·g`)は単項性から**全射**、可逆性から**単射**が自動でつく。
★★★これを知らないと「生成元が非零因子」を手で示す羽目になっていた。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `invertible_principal_equiv` | ★★★★可逆 + 単項 ⟹ `R ≃ₗ I` |
| `principal_of_equiv` | ★★`R ≃ₗ I` ⟹ 単項 |
-/

universe u

namespace ABC3.Found.Arakelov

variable {R : Type u} [CommRing R] (I : Ideal R)

/-- ★★★★**可逆な単項イデアルは自由である**。 -/
theorem invertible_principal_equiv [Module.Invertible R I] (h : Submodule.IsPrincipal I) :
    Nonempty ((R : Type u) ≃ₗ[R] (I : Type u)) := by
  obtain ⟨g, hg⟩ := h
  have hgI : (g : R) ∈ I := by rw [hg]; exact Submodule.mem_span_singleton_self _
  have hsurj : Function.Surjective
      (LinearMap.toSpanSingleton R (I : Type u) ⟨g, hgI⟩) := by
    intro y
    have hy : (y : R) ∈ Submodule.span R {(g : R)} := by rw [← hg]; exact y.2
    obtain ⟨r, hr⟩ := Submodule.mem_span_singleton.1 hy
    exact ⟨r, Subtype.ext (by simpa [LinearMap.toSpanSingleton_apply] using hr)⟩
  exact ⟨LinearEquiv.ofBijective _ (Module.Invertible.bijective_of_surjective hsurj)⟩

theorem principal_of_equiv (e : (R : Type u) ≃ₗ[R] (I : Type u)) :
    Submodule.IsPrincipal I := by
  refine ⟨⟨((e 1 : (I : Type u)) : R), le_antisymm (fun y hy => ?_) ?_⟩⟩
  · have : (⟨y, hy⟩ : (I : Type u)) = e (e.symm ⟨y, hy⟩) := (e.apply_symm_apply _).symm
    refine Submodule.mem_span_singleton.2 ⟨e.symm ⟨y, hy⟩, ?_⟩
    have h2 : e (e.symm ⟨y, hy⟩) = (e.symm ⟨y, hy⟩) • e 1 := by
      rw [← e.map_smul, smul_eq_mul, mul_one]
    exact congrArg Subtype.val (this.trans h2).symm
  · rw [Submodule.span_le, Set.singleton_subset_iff]
    exact (e 1).2


/-! ## ★出典の紐付け(`.src`) -/

def invertible_principal_equiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——可逆な単項イデアルは自由)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
