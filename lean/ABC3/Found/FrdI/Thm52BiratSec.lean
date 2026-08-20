/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm52Rem272

/-!
# [FrdI] Theorem 5.2, (iv) の第 3 段 —— `𝒞^birat` の Frobenius-section と `w_f`

原文 (FrdI p.101):
> for the category C whose objects are objects of C equipped with an F P-path, and

★★本ファイルの中身は 2 つ:

1. **`𝒞^birat` の Frobenius-section** —— `𝒞` の base-Frobenius pair `(𝒫, F)` の像。
   `BaseSection.frobSection`(unit-trivial 型を要求)は `𝒞^birat` には**使えない**ので、
   像で作る。これで `Remark 2.7.2`(等長版)を `𝒞^birat` に当てられる。

2. **`w_f`** —— `F-𝒫-path` つきの対象の間の射 `f : X → Y` に対し
   ```
   w_f := π_X⁻¹ ≫ f^birat ≫ π_Y : ref_X^birat ⟶ ref_Y^birat
   ```
   ★`π` は path が定める `𝒞^birat` の同型(`FPPath.biratIso`)。
   ★`w` は**関手的**(`pathW_comp`, `pathW_id`)である —— 中間の `π_Y⁻¹ ≫ π_Y` が消えるため。

★★★これが `𝒞̃ ⥤ 𝒞^model` の射の単元成分 `u_f` の供給源である:
`w_f` を `Remark 2.7.2` で 3 分解して `w_f = F_n ≫ β_f ≫ α_f` と書き、
`u_f := κ(β_f)` と置く。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

/-! ## ★1. `w_f` —— path が定める参照対象の間の射 -/

/-- ★`FPPath.biratIso` を `BiratCat` の射として名前づけしたもの。

★`HomBirat` のままだと `inv` の型合わせで `instances` 透明度の壁に当たるので、
**射の型を明示した別名**を置く。 -/
noncomputable def FPPath.piHom (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {S : BaseSection P} {X : C} (p : FPPath S X) :
    (toBiratCat P G).obj X ⟶ (toBiratCat P G).obj p.ref :=
  p.biratIso G hiso

instance FPPath.isIso_piHom (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {S : BaseSection P} {X : C} (p : FPPath S X) : IsIso (p.piHom G hiso) :=
  p.biratIso_isIso G hiso

/-- ★★**`w_f := π_X⁻¹ ≫ f^birat ≫ π_Y`** —— path つきの射を、
`𝒞^birat` の中で**参照対象(`𝒫` の対象)の間の射**に移したもの。 -/
noncomputable def pathW (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {S : BaseSection P} {X Y : PathCat S} (f : X ⟶ Y) :
    (toBiratCat P G).obj X.path.ref ⟶ (toBiratCat P G).obj Y.path.ref :=
  inv (X.path.piHom G hiso) ≫ (toBiratCat P G).map f ≫ Y.path.piHom G hiso

/-- ★★★**`w` は合成を保つ** —— 中間の `π_Y ≫ π_Y⁻¹` が消える。 -/
theorem pathW_comp (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {S : BaseSection P} {X Y Z : PathCat S} (f : X ⟶ Y) (g : Y ⟶ Z) :
    pathW G hiso (show X.obj ⟶ Z.obj from f ≫ g) = pathW G hiso f ≫ pathW G hiso g := by
  show inv _ ≫ (toBiratCat P G).map (f ≫ g) ≫ _ = _
  rw [pathW, pathW, (toBiratCat P G).map_comp]
  simp only [Category.assoc, IsIso.hom_inv_id_assoc]

/-- ★`w` は恒等射を保つ。 -/
theorem pathW_id (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {S : BaseSection P} (X : PathCat S) :
    pathW G hiso (𝟙 X) = 𝟙 _ := by
  show inv _ ≫ (toBiratCat P G).map (𝟙 X.obj) ≫ _ = _
  rw [(toBiratCat P G).map_id, Category.id_comp, IsIso.inv_hom_id]

/-! ## ★1b. `w_f` の不変量 -/

@[simp] theorem FPPath.biratDivGp_piHom (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {S : BaseSection P} {X : C} (p : FPPath S X) :
    biratDivGp (p.piHom G hiso) = -p.cls := p.biratDivGp_biratIso G hiso

@[simp] theorem FPPath.biratDeg_piHom (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {S : BaseSection P} {X : C} (p : FPPath S X) :
    biratDeg (p.piHom G hiso) = 1 := p.biratDeg_biratIso G hiso

@[simp] theorem biratDivGp_toBiratMap (G : Frobenioid P) {A B : C} (φ : A ⟶ B) :
    biratDivGp ((toBiratCat P G).map φ) = toGp _ (P.Div φ) :=
  biratDivGp_toHomBirat (P := P) (G := G) φ

@[simp] theorem biratBase_toBiratMap (G : Frobenioid P) {A B : C} (φ : A ⟶ B) :
    biratBase ((toBiratCat P G).map φ) = P.Base φ :=
  biratBase_toHomBirat (P := P) (G := G) φ

@[simp] theorem biratDeg_toBiratMap (G : Frobenioid P) {A B : C} (φ : A ⟶ B) :
    biratDeg ((toBiratCat P G).map φ) = P.degFr φ :=
  biratDeg_toHomBirat (P := P) (G := G) φ

/-- ★★`π_X ≫ w_f = f^birat ≫ π_Y` —— `w_f` の定義そのもの。 -/
theorem pathW_spec (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {S : BaseSection P} {X Y : PathCat S} (f : X ⟶ Y) :
    X.path.piHom G hiso ≫ pathW G hiso f
      = (toBiratCat P G).map f ≫ Y.path.piHom G hiso := by
  rw [pathW, ← Category.assoc, IsIso.hom_inv_id, Category.id_comp]

/-- ★`w_f` の Frobenius 次数は `f` のそれ。 -/
theorem pathW_deg (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {S : BaseSection P} {X Y : PathCat S} (f : X ⟶ Y) :
    biratDeg (pathW G hiso f) = P.degFr f := by
  have h := congrArg (biratPre P G).degFr (pathW_spec G hiso f)
  rw [(biratPre P G).degFr_comp, (biratPre P G).degFr_comp] at h
  have h1 : biratDeg (X.path.piHom G hiso) = 1 := X.path.biratDeg_biratIso G hiso
  have h2 : biratDeg (Y.path.piHom G hiso) = 1 := Y.path.biratDeg_biratIso G hiso
  have h3 : biratDeg ((toBiratCat P G).map f) = P.degFr f :=
    biratDeg_toHomBirat (P := P) (G := G) f
  simp only [biratPre_degFr, h1, h2, h3, mul_one, one_mul] at h
  exact h

/-- ★`w_f` の底は `π` で共役をとったもの。 -/
theorem pathW_base (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {S : BaseSection P} {X Y : PathCat S} (f : X ⟶ Y) :
    biratBase (X.path.piHom G hiso) ≫ biratBase (pathW G hiso f)
      = P.Base f ≫ biratBase (Y.path.piHom G hiso) := by
  have h := congrArg (biratPre P G).Base (pathW_spec G hiso f)
  rw [(biratPre P G).Base_comp, (biratPre P G).Base_comp] at h
  have h3 : biratBase ((toBiratCat P G).map f) = P.Base f :=
    biratBase_toHomBirat (P := P) (G := G) f
  simp only [biratPre_Base, h3] at h
  exact h

/-- ★★★**`w_f` の零因子** —— これが model Frobenioid の条件式そのものになる。

★`Div^gp(π_X) = −cls X`(`FPPath.biratDivGp_biratIso`)を代入すると
```
Φ(Base π_X)(Div^gp w_f) = deg(f)·cls X + Div(f) − Φ(Base f)(cls Y)
```
となり、`Def 5.2, (i)` の関係式の右辺の `Div_B(u_φ)` にちょうど一致する。 -/
theorem pathW_divGp (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {S : BaseSection P} {X Y : PathCat S} (f : X ⟶ Y) :
    Φ.gpMapOn (biratBase (X.path.piHom G hiso)) (biratDivGp (pathW G hiso f))
      = ((P.degFr f : ℕ+) : ℕ) • X.path.cls + toGp _ (P.Div f)
        - Φ.gpMapOn (P.Base f) Y.path.cls := by
  have h := congrArg biratDivGp (pathW_spec G hiso f)
  rw [biratDivGp_comp', biratDivGp_comp'] at h
  have h1 := eq_sub_of_add_eq h
  simp only [FPPath.biratDivGp_piHom, biratDivGp_toBiratMap, biratBase_toBiratMap,
    pathW_deg] at h1
  have e1 : (Φ.gpMapOn (P.Base f) (-Y.path.cls)
        + ((biratDeg (Y.path.piHom G hiso) : ℕ+) : ℕ) • toGp _ (P.Div f))
      - ((P.degFr f : ℕ+) : ℕ) • (-X.path.cls)
      = ((P.degFr f : ℕ+) : ℕ) • X.path.cls + toGp _ (P.Div f)
        - Φ.gpMapOn (P.Base f) Y.path.cls := by
    rw [Y.path.biratDeg_piHom G hiso, map_neg, smul_neg, PNat.one_coe, one_smul]
    abel
  exact h1.trans e1

/-! ## ★2. `𝒞^birat` の Frobenius-section -/

variable [IsConnected D]

/-- ★★★**`𝒞` の Frobenius-section の `𝒞^birat` での像**。

★`BaseSection.frobSection` は unit-trivial 型を要求するので `𝒞^birat` には使えない。
★しかし `𝒞 ⥤ 𝒞^birat` は関手なので、`𝒫` の Frobenius-section をそのまま押し出せる。 -/
noncomputable def biratFrobSection (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    (S : BaseSection P) (Fs : ℕ+ →* SectionEnd S) :
    ℕ+ →* SectionEnd (biratBaseSection G hiso S) where
  toFun n :=
    { app := fun A => ((toBiratCat P G).map
        (((Fs n).app ⟨biratDown P G A.1, A.2⟩ : End _) : _ ⟶ _) : End A.1)
      naturality := by
        rintro A B ⟨f, ⟨f₀, hf, rfl⟩⟩
        show (toBiratCat P G).map f₀ ≫ (toBiratCat P G).map _
          = (toBiratCat P G).map _ ≫ (toBiratCat P G).map f₀
        rw [← (toBiratCat P G).map_comp, ← (toBiratCat P G).map_comp]
        exact congrArg (toBiratCat P G).map
          ((Fs n).naturality (A := (⟨_, A.2⟩ : S.Obj)) (B := (⟨_, B.2⟩ : S.Obj)) ⟨f₀, hf⟩) }
  map_one' := by
    refine SectionEnd.ext ?_
    funext A
    show (toBiratCat P G).map (((Fs 1).app ⟨_, A.2⟩ : End _) : _ ⟶ _) = 𝟙 _
    rw [map_one]
    exact (toBiratCat P G).map_id _
  map_mul' m n := by
    refine SectionEnd.ext ?_
    funext A
    show (toBiratCat P G).map (((Fs (m * n)).app ⟨_, A.2⟩ : End _) : _ ⟶ _)
      = (toBiratCat P G).map (((Fs n).app ⟨_, A.2⟩ : End _) : _ ⟶ _)
        ≫ (toBiratCat P G).map (((Fs m).app ⟨_, A.2⟩ : End _) : _ ⟶ _)
    rw [← (toBiratCat P G).map_comp, map_mul]
    rfl

/-- ★★★**像は再び Frobenius-section** —— 次数・底・Frobenius 型のすべてが移る。 -/
theorem biratFrobSection_isFrobeniusSection (G : Frobenioid P)
    (hiso : ∀ Y : C, IsIsotropic P Y) (S : BaseSection P) {Fs : ℕ+ →* SectionEnd S}
    (hFs : IsFrobeniusSection S Fs) :
    IsFrobeniusSection (biratBaseSection G hiso S) (biratFrobSection G hiso S Fs) where
  degSection n := by
    haveI := (biratBaseSection G hiso S).isConnected_obj
    rw [SectionEnd.deg_eq _ (Classical.arbitrary (biratBaseSection G hiso S).Obj)]
    show biratDeg ((toBiratCat P G).map _) = n
    exact (biratDeg_toHomBirat (G := G) _).trans (hFs.degSection n)
  baseIdentity n A := by
    have hid : biratBase (𝟙 A.1) = P.Base (𝟙 (biratDown P G A.1)) := by
      rw [show (𝟙 A.1) = toHomBirat (𝟙 (biratDown P G A.1)) from
        ((toBiratCat P G).map_id _).symm, biratBase_toHomBirat]
    show biratBase ((toBiratCat P G).map _) = biratBase (𝟙 A.1)
    exact ((biratBase_toHomBirat (G := G)
      (((Fs n).app ⟨biratDown P G A.1, A.2⟩ : End _) : _ ⟶ _)).trans
      (hFs.baseIdentity n ⟨_, A.2⟩)).trans hid.symm
  frobType n A := by
    refine ⟨⟨prop_1_4_i (biratPre P G) _ (fun Y _ => birat_isOfIsotropicType hiso Y),
      birat_isIsometric _⟩, ?_⟩
    show IsIso (biratBase ((toBiratCat P G).map
      (((Fs n).app ⟨biratDown P G A.1, A.2⟩ : End _) : _ ⟶ _)))
    exact (biratBase_toHomBirat (G := G)
      (((Fs n).app ⟨biratDown P G A.1, A.2⟩ : End _) : _ ⟶ _)).symm ▸
      (hFs.frobType n ⟨_, A.2⟩).2

/-! ### ★出典の紐付け(`.src`) -/

/-- ★locator —— `𝒞^birat` の Frobenius-section と `w_f`。 -/
def biratFrobSection.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 101,
    item := "Theorem 5.2, (iv) — 𝒞^birat の Frobenius-section と w_f",
    sectionId := "frdi-thm-5-2" }

end ABC3.Found.FrdI
