/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop44Gl
import ABC3.Found.FrdI.Def27

/-!
# Prop48Sec —— `[FrdI] Proposition 4.8` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.FrdI

universe v u w u2 v2
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

/-! ## ★★★Frobenius-section の移送

原文 (FrdI p.88):
> (iv) If C is of isotropic and pre-model type, then so is Cbirat.

★★`BaseFrobeniusPair` = `sec : BaseSection P` ＋ `frob : ℕ+ →* SectionEnd sec`
＋ `isFrobSection`。★`sec` は上で移したので、残るは `frob` の移送。

★`(S.toBirat Fc).Obj` は `S.Obj` と**台が同じ**なので、
各成分を `toBiratCat` で送るだけ。 -/

variable {G}

/-- ★★**`SectionEnd` の移送** —— 各成分を `toBiratCat` で送る。 -/
noncomputable def SectionEnd.toBirat (Fc : FrobenioidCore (biratPre P G))
    {S : BaseSection P} (ε : SectionEnd S) : SectionEnd (S.toBirat G Fc) where
  app A := (toBiratCat P G).map (ε.app ⟨biratDown P G A.1, A.2⟩)
  naturality {A B} f := by
    obtain ⟨f₀, hf₀, hfe⟩ := f.2
    let A' : S.Obj := ⟨biratDown P G A.1, A.2⟩
    let B' : S.Obj := ⟨biratDown P G B.1, B.2⟩
    have hnat : f₀ ≫ ε.app B' = ε.app A' ≫ f₀ :=
      ε.naturality (⟨f₀, hf₀⟩ : A' ⟶ B')
    show f.1 ≫ (toBiratCat P G).map (ε.app B')
      = (toBiratCat P G).map (ε.app A') ≫ f.1
    exact (congrArg (fun t => t ≫ (toBiratCat P G).map (ε.app B')) hfe.symm).trans
      (((toBiratCat P G).map_comp f₀ (ε.app B')).symm.trans
        ((congrArg (toBiratCat P G).map hnat).trans
          (((toBiratCat P G).map_comp (ε.app A') f₀).trans
            (congrArg (fun t => (toBiratCat P G).map (ε.app A') ≫ t) hfe))))

/-- ★★★**`ℕ+ →* SectionEnd` の移送**。 -/
noncomputable def frobSectionToBirat (Fc : FrobenioidCore (biratPre P G))
    {S : BaseSection P} (Fs : ℕ+ →* SectionEnd S) :
    ℕ+ →* SectionEnd (S.toBirat G Fc) where
  toFun n := SectionEnd.toBirat Fc (Fs n)
  map_one' := by
    refine SectionEnd.ext ?_
    funext A
    show (toBiratCat P G).map ((Fs 1).app _) = (1 : SectionEnd (S.toBirat G Fc)).app A
    rw [Fs.map_one]
    exact (toBiratCat P G).map_id _
  map_mul' m n := by
    refine SectionEnd.ext ?_
    funext A
    show (toBiratCat P G).map ((Fs (m * n)).app _)
      = ((SectionEnd.toBirat Fc (Fs m)) * (SectionEnd.toBirat Fc (Fs n))).app A
    rw [Fs.map_mul]
    exact (toBiratCat P G).map_comp _ _

def frobSectionToBirat.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 88,
    item := "Proposition 4.8, (iv) — Frobenius-section を 𝒞^birat へ移す",
    sectionId := "frdi-prop-4-8" }

/-! ## ★★★★★`Proposition 4.8, (iv)` —— pre-model 型は birat で保たれる -/

/-- ★★★★**Frobenius-section の 3 条も移る**。

★`degSection` は `biratDeg_toHomBirat`、`baseIdentity` は `biratBase_toHomBirat`、
`frobType` は `birat_isFrobeniusType_iff`(すべて在庫)。 -/
theorem isFrobeniusSection_toBirat [IsConnected D] (Fc : FrobenioidCore (biratPre P G))
    {S : BaseSection P} {Fs : ℕ+ →* SectionEnd S} (h : IsFrobeniusSection S Fs) :
    IsFrobeniusSection (S.toBirat G Fc) (frobSectionToBirat Fc Fs) where
  degSection n := by
    haveI := S.isConnected_obj
    let A₀ : S.Obj := Classical.arbitrary S.Obj
    let A : (S.toBirat G Fc).Obj := ⟨A₀.1, A₀.2⟩
    refine (SectionEnd.deg_eq (frobSectionToBirat Fc Fs n) A).trans ?_
    show (biratPre P G).degFr ((toBiratCat P G).map ((Fs n).app _)) = n
    exact (biratDeg_toHomBirat (P := P) (G := G) _).trans
      ((SectionEnd.deg_eq (Fs n) _).symm.trans (h.degSection n))
  baseIdentity n A := by
    show (biratPre P G).Base ((toBiratCat P G).map ((Fs n).app _))
      = (biratPre P G).Base (𝟙 _)
    exact ((biratBase_toHomBirat (P := P) (G := G) _).trans
      ((h.baseIdentity n _).trans (P.Base_id _))).trans ((biratPre P G).Base_id _).symm
  frobType n A :=
    (birat_isFrobeniusType_iff P G _).mpr ⟨(h.frobType n _).1.1, (h.frobType n _).2⟩

/-- ★★★★★**`BaseFrobeniusPair` の移送**。 -/
noncomputable def BaseFrobeniusPair.toBirat [IsConnected D]
    (Fc : FrobenioidCore (biratPre P G)) (BF : BaseFrobeniusPair P) :
    BaseFrobeniusPair (biratPre P G) where
  sec := BF.sec.toBirat G Fc
  frob := frobSectionToBirat Fc BF.frob
  isFrobSection := isFrobeniusSection_toBirat Fc BF.isFrobSection

/-- ★★★★★**[FrdI] Proposition 4.8, (iv)** —— **pre-model 型は birat で保たれる**。

原文 (FrdI p.88):
> (iv) If C is of isotropic and pre-model type, then so is Cbirat.

★isotropic の側は `prop_4_8_i`(在庫)。本定理は pre-model の側。
★★逸脱の開示(**分類 B**): `Fc : FrobenioidCore (biratPre P G)` を要求する
(`BaseSection` の `isPullBack` の条で使う)。
`Proposition 4.4` と同じ birat-Frobenius-normalized 型の仮定である。 -/
theorem isPreModelType_birat [IsConnected D] (Fc : FrobenioidCore (biratPre P G))
    (h : IsPreModelType P) : IsPreModelType (biratPre P G) :=
  h.elim fun BF => ⟨BaseFrobeniusPair.toBirat Fc BF⟩

def prop_4_8_iv.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 88, item := "Proposition 4.8, (iv)",
    sectionId := "frdi-prop-4-8" }

end ABC3.Found.FrdI
