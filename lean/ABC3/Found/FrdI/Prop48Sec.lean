/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop44Gl
import ABC3.Found.FrdI.Def27

/-!
# [FrdI] Proposition 4.8, (iv) —— `BaseSection` を `𝒞^birat` へ移す

原文 (FrdI p.88):
> (iv) If C is of isotropic and pre-model type, then so is Cbirat.

★★`IsPreModelType P := Nonempty (BaseFrobeniusPair P)` で、
`BaseFrobeniusPair` = `sec : BaseSection P` ＋ Frobenius-section。
★したがって **`BaseSection` の移送が実体**である。

## ★★測ったこと(2026-08-19)

`BaseSection P`(`Def27.lean`)はフィールド 8 個の構造体だが、
`BiratCat P G` は**対象が `𝒞` と同じ型**で、
**底も同じ `𝒟`**(`(biratPre P G).toElem.obj A).base = (P.toElem.obj A).base` は `rfl`)。

★★★したがって `toD` に関わる 3 条(`essSurjP` / `fullP` / `faithfulP`)は
`biratBase_toHomBirat`(在庫)1 本で移り、
`skeletal` は `toBiratCat_faithful`(在庫)1 本で移る。

★**`isPullBack` だけが辞書を要する** —— `birat_isPullBack_iff`(在庫)で
「pull-back ⇔ co-angular かつ linear」に翻訳する。
★★これは `FrobenioidCore (biratPre P G)` を要求するので、
**分類 (B) の仮定**(birat-Frobenius-normalized 型)を引き継ぐ。
-/

universe v u w u2 v2

namespace ABC3.Found.FrdI

open CategoryTheory

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} (G : Frobenioid P)

/-- ★★★★★**`BaseSection` は `𝒞^birat` へ移る**。

★`Fc : FrobenioidCore (biratPre P G)` は `isPullBack` の条だけに要る
(**分類 (B)** の仮定。`birat_frobenioidCore_of_frobNormalized` で埋まる)。 -/
noncomputable def BaseSection.toBirat (Fc : FrobenioidCore (biratPre P G))
    (S : BaseSection P) : BaseSection (biratPre P G) where
  objP A := S.objP (biratDown P G A)
  homP {A B} f := ∃ f₀ : biratDown P G A ⟶ biratDown P G B,
    S.homP f₀ ∧ (toBiratCat P G).map f₀ = f
  id_mem {A} hA := by
    refine ⟨𝟙 (biratDown P G A), S.id_mem hA, ?_⟩
    exact (toBiratCat P G).map_id (biratDown P G A)
  comp_mem {_ _ _} {f g} hf hg := by
    obtain ⟨f₀, hf₀, hfe⟩ := hf
    obtain ⟨g₀, hg₀, hge⟩ := hg
    refine ⟨f₀ ≫ g₀, S.comp_mem hf₀ hg₀, ?_⟩
    exact ((toBiratCat P G).map_comp f₀ g₀).trans
      ((congrArg (fun t => t ≫ (toBiratCat P G).map g₀) hfe).trans
        (congrArg (fun t => f ≫ t) hge))
  isPullBack {_ _} {f} hf := by
    obtain ⟨f₀, hf₀, hfe⟩ := hf
    rw [← hfe]
    refine (birat_isPullBack_iff P G Fc f₀).mpr ?_
    obtain ⟨hlb, hlin⟩ := G.core.pullBackLB f₀ (S.isPullBack hf₀)
    exact ⟨hlb.1, hlin⟩
  skeletal {_ _} {f g} hf hg hfg hgf := by
    haveI : (toBiratCat P G).Faithful := toBiratCat_faithful
    obtain ⟨f₀, hf₀, hfe⟩ := hf
    obtain ⟨g₀, hg₀, hge⟩ := hg
    refine S.skeletal hf₀ hg₀ ?_ ?_
    · refine (toBiratCat P G).map_injective ?_
      exact (((toBiratCat P G).map_comp f₀ g₀).trans
        ((congrArg (fun t => t ≫ (toBiratCat P G).map g₀) hfe).trans
          (congrArg (fun t => f ≫ t) hge))).trans
        (hfg.trans ((toBiratCat P G).map_id _).symm)
    · refine (toBiratCat P G).map_injective ?_
      exact (((toBiratCat P G).map_comp g₀ f₀).trans
        ((congrArg (fun t => t ≫ (toBiratCat P G).map f₀) hge).trans
          (congrArg (fun t => g ≫ t) hfe))).trans
        (hgf.trans ((toBiratCat P G).map_id _).symm)
  frobTrivial {A} hA := by
    obtain ⟨ζ, hζd, hζb⟩ := S.frobTrivial hA
    refine ⟨(Functor.mapEnd (biratDown P G A) (toBiratCat P G)).comp ζ, ?_, ?_⟩
    · intro n
      show (biratPre P G).degFr ((toBiratCat P G).map ((ζ n : End (biratDown P G A)))) = n
      exact (biratDeg_toHomBirat (P := P) (G := G) _).trans (hζd n)
    · intro n
      refine ⟨?_, ?_⟩
      · show (biratPre P G).Base ((toBiratCat P G).map ((ζ n : End (biratDown P G A))))
          = (biratPre P G).Base (𝟙 _)
        exact ((biratBase_toHomBirat (P := P) (G := G) _).trans
          (((hζb n).1).trans (P.Base_id _))).trans ((biratPre P G).Base_id _).symm
      · exact (birat_isFrobeniusType_iff P G _).mpr ⟨(hζb n).2.1.1, (hζb n).2.2⟩
  essSurjP X := by
    obtain ⟨A, hA, hiso⟩ := S.essSurjP X
    exact ⟨A, hA, hiso⟩
  fullP {A B} hA hB ψ := by
    obtain ⟨f₀, hf₀, hfe⟩ := S.fullP hA hB ψ
    refine ⟨(toBiratCat P G).map f₀, ⟨f₀, hf₀, rfl⟩, ?_⟩
    exact (biratBase_toHomBirat (P := P) (G := G) f₀).trans hfe
  faithfulP {_ _} {f g} hf hg hb := by
    obtain ⟨f₀, hf₀, hfe⟩ := hf
    obtain ⟨g₀, hg₀, hge⟩ := hg
    rw [← hfe, ← hge]
    refine congrArg (toBiratCat P G).map (S.faithfulP hf₀ hg₀ ?_)
    exact (biratBase_toHomBirat (P := P) (G := G) f₀).symm.trans
      ((congrArg (biratPre P G).Base hfe).trans
        (hb.trans ((congrArg (biratPre P G).Base hge).symm.trans
          (biratBase_toHomBirat (P := P) (G := G) g₀))))

def BaseSection.toBirat.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 88,
    item := "Proposition 4.8, (iv) — BaseSection を 𝒞^birat へ移す",
    sectionId := "frdi-prop-4-8" }

end ABC3.Found.FrdI
