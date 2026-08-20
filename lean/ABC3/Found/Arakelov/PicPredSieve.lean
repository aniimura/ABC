import ABC3.Found.Arakelov.PicHcompatImg

/-!
# Arakelov (B2) 第 236 ブロック —— **「アフィン」を任意の述語に緩める**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★第 189 の仮定は**強すぎた**

第 189 ブロック `locallyBijective_of_bijective_on_affine` は
「**すべての**アフィン開で全単射」を要求する。

★ところが引き戻しの比較射が全単射になるのは

    V はアフィン ∧ V ≤ f⁻¹ᵁ B(B は `Y` のアフィン開)∧ V で両辺が自明

を満たす `V` **だけ**である。`X` の勝手なアフィン開 `A` に対して
`f(A)` を含む `Y` のアフィン開が在るとは限らない。

★★逃げ道は第 234 と同じ——**述語をパラメータにする**。
必要なのは「点ごとに細かく取れる」ことだけである:

    ∀ U, ∀ x ∈ U, ∃ A, P A ∧ x ∈ A ∧ A ≤ U

★★★証明は第 189 と 1 文字も変わらない。

| 定義・定理 | 内容 |
|---|---|
| `predSieve` | ★`P` を満たす開集合で生成される篩 |
| `predSieve_mem` | ★★点ごとに取れれば被覆篩 |
| `locallyBijective_of_bijective_on_pred` | ★★★★**局所全単射の判定(一般形)** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X : Scheme.{u}}

/-- ★**述語 `P` を満たす開集合で生成される篩**。 -/
def predSieve (P : X.Opens → Prop) (U : X.Opens) : Sieve U where
  arrows V _ := ∃ A : X.Opens, P A ∧ V ≤ A ∧ A ≤ U
  downward_closed := by
    rintro V W i ⟨A, hPA, hVA, hAU⟩ j
    exact ⟨A, hPA, le_trans (leOfHom j) hVA, hAU⟩

/-- ★★**点ごとに取れれば被覆篩である**。 -/
theorem predSieve_mem (P : X.Opens → Prop)
    (hP : ∀ (U : X.Opens) (x : X), x ∈ U → ∃ A : X.Opens, P A ∧ x ∈ A ∧ A ≤ U)
    (U : X.Opens) : predSieve P U ∈ (Opens.grothendieckTopology X) U := by
  intro x hx
  obtain ⟨A, hPA, hxA, hAU⟩ := hP U x hx
  exact ⟨A, homOfLE hAU, ⟨A, hPA, le_rfl, hAU⟩, hxA⟩

/-- ★★★★**局所全単射の判定(一般形)**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★第 189 の「アフィン開」を任意の述語 `P` に緩めたもの。 -/
theorem locallyBijective_of_bijective_on_pred
    {F G : (X.Opens)ᵒᵖ ⥤ AddCommGrpCat.{u}} (f : F ⟶ G) (P : X.Opens → Prop)
    (hP : ∀ (U : X.Opens) (x : X), x ∈ U → ∃ A : X.Opens, P A ∧ x ∈ A ∧ A ≤ U)
    (h : ∀ A : X.Opens, P A → Function.Bijective (f.app (op A))) :
    Presheaf.IsLocallyInjective (Opens.grothendieckTopology X) f ∧
      Presheaf.IsLocallySurjective (Opens.grothendieckTopology X) f := by
  constructor
  · refine isLocallyInjective_of_coverSieve _ f (fun U x y hxy => ?_)
    refine ⟨predSieve P U, predSieve_mem P hP U, ?_⟩
    rintro V i ⟨A, hPA, hVA, hAU⟩
    have hA : F.map (homOfLE hAU).op x = F.map (homOfLE hAU).op y := by
      refine (h A hPA).1 ?_
      rw [NatTrans.naturality_apply f (homOfLE hAU).op x,
        NatTrans.naturality_apply f (homOfLE hAU).op y, hxy]
    have h2 := congrArg (fun z => F.map (homOfLE hVA).op z) hA
    simp only [presheafRes_trans] at h2
    exact h2
  · refine isLocallySurjective_of_cover _ f (fun U s => ?_)
    refine ⟨predSieve P U, predSieve_mem P hP U, ?_⟩
    rintro V i ⟨A, hPA, hVA, hAU⟩
    obtain ⟨tA, htA⟩ := (h A hPA).2 (G.map (homOfLE hAU).op s)
    refine ⟨F.map (homOfLE hVA).op tA, ?_⟩
    rw [NatTrans.naturality_apply f (homOfLE hVA).op tA, htA, presheafRes_trans]
    exact rfl

/-! ## ★出典の紐付け(`.src`) -/

def locallyBijective_of_bijective_on_pred.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——局所全単射の判定を任意の述語で)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
