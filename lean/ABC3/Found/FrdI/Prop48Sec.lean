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
  id_mem {A} hA := ⟨𝟙 _, S.id_mem hA, (toBiratCat P G).map_id _⟩
  comp_mem {_ _ _} {f g} hf hg := by
    obtain ⟨f₀, hf₀, hfe⟩ := hf
    obtain ⟨g₀, hg₀, hge⟩ := hg
    refine ⟨f₀ ≫ g₀, S.comp_mem hf₀ hg₀, ?_⟩
    rw [(toBiratCat P G).map_comp, hfe, hge]
  isPullBack {_ _} {f} hf := by
    obtain ⟨f₀, hf₀, hfe⟩ := hf
    rw [← hfe]
    refine (birat_isPullBack_iff G Fc f₀).mpr ?_
    obtain ⟨hlb, hlin⟩ := G.core.pullBackLB f₀ (S.isPullBack hf₀)
    exact ⟨hlb.1, hlin⟩
  skeletal {_ _} {f g} hf hg hfg hgf := by
    haveI : (toBiratCat P G).Faithful := toBiratCat_faithful
    obtain ⟨f₀, hf₀, hfe⟩ := hf
    obtain ⟨g₀, hg₀, hge⟩ := hg
    refine S.skeletal hf₀ hg₀ ?_ ?_
    · refine (toBiratCat P G).map_injective ?_
      rw [(toBiratCat P G).map_comp, hfe, hge, hfg, (toBiratCat P G).map_id]
    · refine (toBiratCat P G).map_injective ?_
      rw [(toBiratCat P G).map_comp, hge, hfe, hgf, (toBiratCat P G).map_id]
  frobTrivial {A} hA := by
    obtain ⟨ζ, hζd, hζb⟩ := S.frobTrivial hA
    refine ⟨(Functor.mapEnd (biratDown P G A) (toBiratCat P G)).comp ζ, ?_, ?_⟩
    · intro n
      show (biratPre P G).degFr ((toBiratCat P G).map ((ζ n : End (biratDown P G A)))) = n
      rw [show (biratPre P G).degFr ((toBiratCat P G).map ((ζ n : End (biratDown P G A))))
        = P.degFr ((ζ n : End (biratDown P G A))) from biratDeg_toHomBirat _]
      exact hζd n
    · intro n
      refine ⟨?_, ?_⟩
      · show (biratPre P G).Base ((toBiratCat P G).map ((ζ n : End (biratDown P G A))))
          = (biratPre P G).Base (𝟙 _)
        rw [(biratPre P G).Base_id,
          show (biratPre P G).Base ((toBiratCat P G).map ((ζ n : End (biratDown P G A))))
            = P.Base ((ζ n : End (biratDown P G A))) from biratBase_toHomBirat _,
          show P.Base ((ζ n : End (biratDown P G A))) = P.Base (𝟙 _) from (hζb n).1,
          P.Base_id]
      · exact (birat_isFrobeniusType_iff G _).mpr ⟨(hζb n).2.1.1, (hζb n).2.2⟩
  essSurjP X := by
    obtain ⟨A, hA, hiso⟩ := S.essSurjP X
    exact ⟨A, hA, hiso⟩
  fullP {A B} hA hB ψ := by
    obtain ⟨f₀, hf₀, hfe⟩ := S.fullP hA hB ψ
    refine ⟨(toBiratCat P G).map f₀, ⟨f₀, hf₀, rfl⟩, ?_⟩
    rw [show (biratPre P G).Base ((toBiratCat P G).map f₀) = P.Base f₀ from
      biratBase_toHomBirat f₀]
    exact hfe
  faithfulP {_ _} {f g} hf hg hb := by
    obtain ⟨f₀, hf₀, hfe⟩ := hf
    obtain ⟨g₀, hg₀, hge⟩ := hg
    rw [← hfe, ← hge]
    refine congrArg (toBiratCat P G).map (S.faithfulP hf₀ hg₀ ?_)
    rw [← biratBase_toHomBirat (P := P) (G := G) f₀,
      ← biratBase_toHomBirat (P := P) (G := G) g₀, hfe, hge]
    exact hb

def BaseSection.toBirat.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 88,
    item := "Proposition 4.8, (iv) — BaseSection を 𝒞^birat へ移す",
    sectionId := "frdi-prop-4-8" }

end ABC3.Found.FrdI
