import ABC3.Found.FrdI.Prop15
import ABC3.Found.FrdI.Witness

/-!
# 捻れ積 `𝔽_Φ ⋉ G` —— `[FrdI] Proposition 1.14, (iii)` の反例の候補

★`Gap/FrdI/Section1.lean` の `Gap_1_14_iii` を ③(`sourceGap`)に上げるための構成。

## ★なぜ単なる直積では駄目か

合成を `(f, g) ≫ (f′, g′) = (f ≫ f′, g · g′)` と定めると、
`Frobenius-normalized` の等式 `ζ ≫ α^d = α ≫ ζ` の `G` 成分は
`g_ζ · g_α^d = g_α · g_ζ` となり、★**`g_α^{d−1} = 1` をすべての `d` について
強制する**ので `G` が自明になってしまう。

## ★★捻れ積

合成を **Frobenius 次数で捻る**:
```
(f, g) ≫ (f′, g′) = (f ≫ f′, g ^ degFr(f′) · g′)
```
`Base`・`Div`・`degFr` は `G` 成分を無視する(射影関手 `TwProj`)。

★**`Frobenius-normalized` が自動で成り立つ**:
`α = (𝟙, a)`(次数 1)、`ζ = (z, 1)`(次数 `d`)に対し
`ζ ≫ α^d = (z, 1^1 · a^d) = (z, a^d)`、`α ≫ ζ = (z, a^d · 1) = (z, a^d)`。

★★**そして `ζ` は mono でない**: `(f₀, u) ≫ (z, 1) = (f₀ ≫ z, u^d)` なので、
`u^d = v^d` から `u = v` は出ない。★**`G` が `d`-捻れを持てば破れる。**

## ★残る仕事

★**この捻れ積が `Definition 1.3` の全条件を満たすことの検証**である。
`𝔽_Φ` が Frobenioid であることは `Proposition 1.5`(`elemFrob_isFrobenioid`)で
証明済みなので、各条件に `G` 成分の処理を足す形になると見ている(★未検証)。
-/

namespace ABC3.Check.FrdI

open CategoryTheory ABC3.Found.FrdI

universe v u w t

variable {D : Type u} [Category.{v} D] {Φ : MonoidOn.{v, u, w} D}
  {G : Type t} [CommGroup G]

/-- ★**捻れ積の対象** —— `𝔽_Φ` の対象そのもの(型シノニムとして包む)。

★`G` は本体に現れないが、**圏構造が `G` に依存する**ので型に含める。 -/
structure TwObj (Φ : MonoidOn.{v, u, w} D) (G : Type t) [CommGroup G] where
  /-- 元の `𝔽_Φ` の対象 -/
  ofElem : ElemFrobCat Φ

/-- ★**捻れ積の射** —— `𝔽_Φ` の射に `G` の元を添えたもの。 -/
@[ext]
structure TwHom (A B : TwObj Φ G) where
  /-- `𝔽_Φ` 成分 -/
  hom : A.ofElem ⟶ B.ofElem
  /-- `G` 成分 -/
  unit : G

/-- ★★**捻れ積の圏構造** —— 合成の `G` 成分を **Frobenius 次数で捻る**。

★結合律は `(g ≫ h).deg = h.deg * g.deg`(`𝔽_Φ` では `rfl`)と
`G` の可換性から出る。 -/
instance twCategory : Category (TwObj Φ G) where
  Hom A B := TwHom A B
  id A := ⟨𝟙 A.ofElem, 1⟩
  comp f g := ⟨f.hom ≫ g.hom, f.unit ^ ((g.hom.deg : ℕ+) : ℕ) * g.unit⟩
  id_comp f := by
    refine TwHom.ext (Category.id_comp _) ?_
    show (1 : G) ^ ((f.hom.deg : ℕ+) : ℕ) * f.unit = f.unit
    rw [one_pow, one_mul]
  comp_id := fun {_ Y} f => by
    refine TwHom.ext (Category.comp_id _) ?_
    show f.unit ^ ((ElemFrobCat.Hom.deg (𝟙 Y.ofElem) : ℕ+) : ℕ) * (1 : G) = f.unit
    rw [show ElemFrobCat.Hom.deg (𝟙 Y.ofElem) = 1 from rfl]
    simp
  assoc f g h := by
    refine TwHom.ext (Category.assoc _ _ _) ?_
    show (f.unit ^ ((g.hom.deg : ℕ+) : ℕ) * g.unit) ^ ((h.hom.deg : ℕ+) : ℕ) * h.unit
      = f.unit ^ (((g.hom ≫ h.hom).deg : ℕ+) : ℕ)
          * (g.unit ^ ((h.hom.deg : ℕ+) : ℕ) * h.unit)
    rw [show (g.hom ≫ h.hom).deg = h.hom.deg * g.hom.deg from rfl, PNat.mul_coe,
      mul_pow, ← pow_mul, mul_comm ((h.hom.deg : ℕ+) : ℕ), mul_assoc]

@[simp] theorem twComp_hom {A B E : TwObj Φ G} (f : A ⟶ B) (g : B ⟶ E) :
    (f ≫ g).hom = f.hom ≫ g.hom := rfl

@[simp] theorem twComp_unit {A B E : TwObj Φ G} (f : A ⟶ B) (g : B ⟶ E) :
    (f ≫ g).unit = f.unit ^ ((g.hom.deg : ℕ+) : ℕ) * g.unit := rfl

@[simp] theorem twId_hom (A : TwObj Φ G) : (𝟙 A : A ⟶ A).hom = 𝟙 A.ofElem := rfl

@[simp] theorem twId_unit (A : TwObj Φ G) : (𝟙 A : A ⟶ A).unit = 1 := rfl

/-- ★**射影関手 `𝔽_Φ ⋉ G → 𝔽_Φ`** —— `G` 成分を落とす。

★★**これが pre-Frobenioid 構造を与える** —— `Base`・`Div`・`degFr` は
すべて `G` 成分を無視する。 -/
def twProj (Φ : MonoidOn.{v, u, w} D) (G : Type t) [CommGroup G] :
    TwObj Φ G ⥤ ElemFrobCat Φ where
  obj A := A.ofElem
  map f := f.hom

/-! ## ★★核心 —— 次数 `d` の射は `G` の `d`-捻れを潰す -/

/-- ★★**捻れ積では、次数 `d` の射での右合成が `G` 成分を `d` 乗する**。 -/
theorem twComp_unit_pow {A B E : TwObj Φ G} (f : A ⟶ B) (ζ : B ⟶ E)
    (hζ : ζ.unit = 1) : (f ≫ ζ).unit = f.unit ^ ((ζ.hom.deg : ℕ+) : ℕ) := by
  rw [twComp_unit, hζ, mul_one]

/-- ★★★**次数 `d` の射は mono でない** —— `G` に `d`-捻れがあるとき。

★`u ≠ 1` かつ `u ^ d = 1` なる `u ∈ G` があれば、
`(𝟙, 1) ≫ ζ = (𝟙, u) ≫ ζ` だが `(𝟙, 1) ≠ (𝟙, u)`。

★★**これが `Proposition 1.14, (iii)` の `⟸` を破る機構である** ——
そこは「prime-Frobenius 射が FSMI(したがって mono)」を要求する。 -/
theorem not_mono_twist {B E : TwObj Φ G} (ζ : B ⟶ E) (hζ : ζ.unit = 1)
    (u : G) (hu : u ≠ 1) (hup : u ^ ((ζ.hom.deg : ℕ+) : ℕ) = 1) : ¬ Mono ζ := by
  intro hm
  have hsq : (⟨𝟙 B.ofElem, 1⟩ : TwHom B B) ≫ ζ = (⟨𝟙 B.ofElem, u⟩ : TwHom B B) ≫ ζ := by
    refine TwHom.ext rfl ?_
    rw [twComp_unit, twComp_unit, hζ, mul_one, mul_one, one_pow, hup]
  have := (cancel_mono ζ).mp hsq
  exact hu (congrArg TwHom.unit this).symm

/-! ## ★捻れ積の `PreFrobenioid` 構造

★`𝔽_Φ` の側の結果に **`G` 成分の処理を足す**だけで出る。
-/

/-- ★**捻れ積も totally epimorphic**。

★`𝔽_Φ` 成分は `isTotallyEpimorphic_elemFrobCat`、
★**`G` 成分は群の消約**(指数が一致するので左から消せる)。 -/
theorem twTotallyEpimorphic (hD : IsTotallyEpimorphic D)
    (hint : ∀ A : D, IsIntegralMonoid (Φ.val A)) :
    IsTotallyEpimorphic (TwObj Φ G) := by
  intro A B f
  haveI : Epi f.hom := isTotallyEpimorphic_elemFrobCat hD hint _ _ f.hom
  refine ⟨fun {Z} h₁ h₂ hh => ?_⟩
  have hhom : h₁.hom = h₂.hom := (cancel_epi f.hom).mp (congrArg TwHom.hom hh)
  refine TwHom.ext hhom ?_
  have hu := congrArg TwHom.unit hh
  rw [twComp_unit, twComp_unit, hhom] at hu
  exact mul_left_cancel hu

/-- ★**捻れ積も connected** —— `𝔽_Φ` の証明をそのまま写す。 -/
theorem twIsConnected [IsConnected D] : IsConnected (TwObj Φ G) := by
  obtain ⟨d₀⟩ := (inferInstance : Nonempty D)
  refine IsConnected.of_induct (j₀ := (⟨⟨d₀⟩⟩ : TwObj Φ G)) ?_
  intro p hp0 hstep A
  have key : ∀ d : D, (⟨⟨d⟩⟩ : TwObj Φ G) ∈ p :=
    induct_on_objects (J := D) {d | (⟨⟨d⟩⟩ : TwObj Φ G) ∈ p} hp0
      (fun {d₁ d₂} f => hstep
        (⟨(⟨f, 0, 1⟩ : (⟨d₁⟩ : ElemFrobCat Φ) ⟶ (⟨d₂⟩ : ElemFrobCat Φ)), 1⟩ :
          (⟨⟨d₁⟩⟩ : TwObj Φ G) ⟶ ⟨⟨d₂⟩⟩))
  exact key A.ofElem.base

/-- ★★**捻れ積の pre-Frobenioid 構造**。

★`𝔽_Φ` の pre-Frobenioid 構造に射影関手を前合成するだけ ——
★**`Base`・`Div`・`degFr` はすべて `G` 成分を無視する。** -/
def twPreFrobenioid (hD : IsTotallyEpimorphic D) [IsConnected D]
    (hpd : ∀ A : D, IsPreDivisorial (Φ.val A)) :
    PreFrobenioid (TwObj Φ G) Φ.charOn where
  toElem := twProj Φ G ⋙ elemFrobToChar Φ
  divisorial := Φ.charOn_isDivisorialOn hpd
  totEpiC := twTotallyEpimorphic hD (fun A => (hpd A).1)
  totEpiD := hD
  connectedC := twIsConnected
  connectedD := inferInstance

/-! ## ★同型は `𝔽_Φ` 成分だけで決まる

★★**これが「捻れ積も isotropic 型である」ことの土台**である ——
`Div = 0` の pre-step は `𝔽_Φ` で同型になり、`G` 成分は自動で逆が取れる。
-/

/-- ★★**`𝔽_Φ` 成分が同型なら、捻れ積でも同型**。

★逆射の `G` 成分は `f.unit⁻¹`。★**同型の次数は 1**(`ElemFrobCat.isIso_iff`)
なので指数が消え、両側の等式がそろう。 -/
theorem twIsIso_of {A B : TwObj Φ G} (f : A ⟶ B) (h : IsIso f.hom) : IsIso f := by
  haveI := h
  obtain ⟨-, -, hdeg⟩ := (ElemFrobCat.isIso_iff f.hom).mp h
  have hdegi : ElemFrobCat.Hom.deg (inv f.hom) = 1 := by
    have h1 := congrArg ElemFrobCat.Hom.deg (IsIso.inv_hom_id f.hom)
    rw [ElemFrobCat.degFr_comp, hdeg, one_mul] at h1
    exact h1
  refine ⟨⟨⟨inv f.hom, f.unit⁻¹⟩, ?_, ?_⟩⟩
  · refine TwHom.ext (IsIso.hom_inv_id f.hom) ?_
    show f.unit ^ ((ElemFrobCat.Hom.deg (inv f.hom) : ℕ+) : ℕ) * f.unit⁻¹ = (1 : G)
    rw [hdegi]
    simp
  · refine TwHom.ext (IsIso.inv_hom_id f.hom) ?_
    show f.unit⁻¹ ^ ((ElemFrobCat.Hom.deg f.hom : ℕ+) : ℕ) * f.unit = (1 : G)
    rw [hdeg]
    simp

/-! ## ★不変量は `𝔽_Φ` 成分だけで決まる(橋渡し)

★★**これが「`G` 成分を無視する」の正確な意味**である。 -/

section Bridge

variable (hD : IsTotallyEpimorphic D) [IsConnected D]
  (hpd : ∀ A : D, IsPreDivisorial (Φ.val A))

@[simp] theorem twBase {A B : TwObj Φ G} (f : A ⟶ B) :
    (twPreFrobenioid (G := G) hD hpd).Base f = f.hom.base := rfl

@[simp] theorem twDegFr {A B : TwObj Φ G} (f : A ⟶ B) :
    (twPreFrobenioid (G := G) hD hpd).degFr f = f.hom.deg := rfl

@[simp] theorem twDiv {A B : TwObj Φ G} (f : A ⟶ B) :
    (twPreFrobenioid (G := G) hD hpd).Div f = toChar f.hom.div := rfl

/-- ★★**捻れ積も isotropic 型** —— `Div = 0` の pre-step は同型。

★`𝔽_Φ` の側(`elemFrob_isotropic`)で `hom` が同型になり、
`twIsIso_of` が `G` 成分の逆を自動で作る。 -/
theorem twIsotropic (A : TwObj Φ G) : IsIsotropic (twPreFrobenioid (G := G) hD hpd) A := by
  intro Dd φ hisom hstep
  refine twIsIso_of φ ?_
  exact elemFrob_isotropic Φ hD hpd A.ofElem Dd.ofElem φ.hom hisom hstep

/-- ★★**捻れ積では全射が co-angular**(`Proposition 1.4, (i)`)。

★★**これで `Definition 1.3` の (iii)(a)(b) が自明になる。** -/
theorem twCoAngular {A B : TwObj Φ G} (φ : A ⟶ B) :
    IsCoAngular (twPreFrobenioid (G := G) hD hpd) φ :=
  prop_1_4_i _ φ (fun Y _ => twIsotropic hD hpd Y)

/-- ★**捻れ積でも pre-step は mono**(`Definition 1.3, (v), (a)`)。

★`𝔽_Φ` 成分は `elemFrob_preStepMono`、★**`G` 成分は次数 1 なので単純に消約**。 -/
theorem twPreStepMono {A B : TwObj Φ G} (φ : A ⟶ B)
    (hφ : IsPreStep (twPreFrobenioid (G := G) hD hpd) φ) : Mono φ := by
  haveI : Mono φ.hom := elemFrob_preStepMono Φ hD hpd φ.hom hφ
  refine ⟨fun {Z} f g hh => ?_⟩
  have hhom : f.hom = g.hom := (cancel_mono φ.hom).mp (congrArg TwHom.hom hh)
  refine TwHom.ext hhom ?_
  have hu := congrArg TwHom.unit hh
  rw [twComp_unit, twComp_unit, show ElemFrobCat.Hom.deg φ.hom = 1 from hφ.1] at hu
  simpa using mul_right_cancel hu

/-- ★**各次数の Frobenius 型自己射がある**(`Definition 1.3, (ii)` の存在部分)。

★`(⟨𝟙, 0, n⟩, 1)` を取ればよい。★**co-angular は `twCoAngular` が無条件に与える。** -/
theorem twFrobDegSurj (A : TwObj Φ G) (n : ℕ+) :
    ∃ (B : TwObj Φ G) (φ : A ⟶ B),
      IsFrobeniusType (twPreFrobenioid (G := G) hD hpd) φ ∧
      (twPreFrobenioid (G := G) hD hpd).degFr φ = n := by
  refine ⟨A, ⟨⟨𝟙 A.ofElem.base, 0, n⟩, 1⟩, ⟨⟨twCoAngular hD hpd _, ?_⟩, ?_⟩, rfl⟩
  · show toChar (0 : Φ.val A.ofElem.base) = 0
    simp
  · show IsIso (𝟙 A.ofElem.base)
    infer_instance

/-- ★**(iii)(a)** co-angular は合成で閉じる —— ★**全射が co-angular なので自明**。 -/
theorem twCoAngularComp {A B E : TwObj Φ G} (ψ : A ⟶ B) (φ : B ⟶ E)
    (_ : IsCoAngular (twPreFrobenioid (G := G) hD hpd) ψ)
    (_ : IsCoAngular (twPreFrobenioid (G := G) hD hpd) φ) :
    IsCoAngular (twPreFrobenioid (G := G) hD hpd) (ψ ≫ φ) :=
  twCoAngular hD hpd _

/-- ★**(vii)(b)** isotropic は下流へ伝播 —— ★**全対象が isotropic なので自明**。 -/
theorem twIsotropicClosed {A B : TwObj Φ G} (_ : A ⟶ B)
    (_ : IsIsotropic (twPreFrobenioid (G := G) hD hpd) A) :
    IsIsotropic (twPreFrobenioid (G := G) hD hpd) B :=
  twIsotropic hD hpd B

/-- ★**(ii) の本質的一意性** —— 同じ次数の Frobenius 型射は同型を除いて一意。

★`𝔽_Φ` 側の `β₀` に、★**`G` 成分を `φ.unit⁻¹ * ψ.unit` と取れば四角形が閉じる**
(`β₀` は同型なので次数 1)。 -/
theorem twFrobType_hom {A B : TwObj Φ G} (φ : A ⟶ B)
    (hφ : IsFrobeniusType (twPreFrobenioid (G := G) hD hpd) φ) :
    IsFrobeniusType (elemPreFrobenioid Φ hD hpd) φ.hom :=
  ⟨⟨prop_1_4_i _ φ.hom (fun Y _ => elemFrob_isotropic Φ hD hpd Y), hφ.1.2⟩, hφ.2⟩

theorem twFrobDegUniq (A B E : TwObj Φ G) (φ : A ⟶ B) (ψ : A ⟶ E)
    (hφ : IsFrobeniusType (twPreFrobenioid (G := G) hD hpd) φ)
    (hψ : IsFrobeniusType (twPreFrobenioid (G := G) hD hpd) ψ)
    (hd : (twPreFrobenioid (G := G) hD hpd).degFr φ
        = (twPreFrobenioid (G := G) hD hpd).degFr ψ) :
    ∃ β : B ⟶ E, IsIso β ∧ φ ≫ β = ψ := by
  obtain ⟨β₀, hβ₀iso, hβ₀⟩ :=
    elemFrob_frobDegUniq Φ hD hpd A.ofElem B.ofElem E.ofElem φ.hom ψ.hom
      (twFrobType_hom hD hpd φ hφ) (twFrobType_hom hD hpd ψ hψ) hd
  haveI := hβ₀iso
  obtain ⟨-, -, hdβ₀⟩ := (ElemFrobCat.isIso_iff β₀).mp hβ₀iso
  refine ⟨⟨β₀, φ.unit⁻¹ * ψ.unit⟩, twIsIso_of _ hβ₀iso, ?_⟩
  refine TwHom.ext hβ₀ ?_
  show φ.unit ^ ((ElemFrobCat.Hom.deg β₀ : ℕ+) : ℕ) * (φ.unit⁻¹ * ψ.unit) = ψ.unit
  rw [hdβ₀]
  simp

/-- ★`ζ n := (⟨𝟙, 0, n⟩, 1)` —— Frobenius-trivial 性を与えるモノイド準同型。

★★**`End` の乗法は `x * y = y ≫ x`** なので、`map_mul'` は
`ζ (n*m) = ζ m ≫ ζ n` を示すことになる。 -/
def twZeta (Y : D) : ℕ+ →* End (⟨⟨Y⟩⟩ : TwObj Φ G) where
  toFun n := ⟨⟨𝟙 Y, 0, n⟩, 1⟩
  map_one' := TwHom.ext (ElemFrobCat.Hom.ext rfl rfl rfl) rfl
  map_mul' n m := by
    refine TwHom.ext ?_ ?_
    · refine ElemFrobCat.Hom.ext ?_ ?_ rfl
      · show 𝟙 Y = 𝟙 Y ≫ 𝟙 Y
        simp
      · show (0 : Φ.val Y) = Φ.map (𝟙 Y) 0 + ((n : ℕ+) : ℕ) • (0 : Φ.val Y)
        simp
    · show (1 : G) = (1 : G) ^ ((n : ℕ+) : ℕ) * 1
      simp

/-- ★**(i)(a)** 各底に Frobenius-trivial 対象がある。 -/
theorem twBaseSurj (Y : D) :
    ∃ A : TwObj Φ G, IsFrobeniusTrivial (twPreFrobenioid (G := G) hD hpd) A ∧
      Nonempty (((twPreFrobenioid (G := G) hD hpd).toElem.obj A).base ≅ Y) := by
  refine ⟨⟨⟨Y⟩⟩, ⟨twZeta (G := G) Y, fun n => rfl, fun n => ?_⟩, ⟨Iso.refl _⟩⟩
  exact ⟨rfl, ⟨twCoAngular hD hpd _, by show toChar (0 : Φ.val Y) = 0; simp⟩,
    by show IsIso (𝟙 Y); infer_instance⟩

/-- ★**(iii)(b)** —— ★**全射が co-angular なので自明**。 -/
theorem twCoAngularOfPreStep {A' A : TwObj Φ G} (α : A' ⟶ A)
    (_ : IsCoAngular (twPreFrobenioid (G := G) hD hpd) α)
    (_ : IsPreStep (twPreFrobenioid (G := G) hD hpd) α) (φ : A' ⟶ A) :
    IsCoAngular (twPreFrobenioid (G := G) hD hpd) φ :=
  twCoAngular hD hpd φ

/-- ★**(vii)(a)** isotropic hull は `𝟙` —— ★**全対象が isotropic だから**。 -/
theorem twIsotropicHullExists (A : TwObj Φ G) :
    ∃ (B : TwObj Φ G) (φ : A ⟶ B),
      IsIsotropicHull (twPreFrobenioid (G := G) hD hpd) φ :=
  ⟨A, 𝟙 A, (twPreFrobenioid (G := G) hD hpd).Div_id A,
    isPreStep_id _ A, twIsotropic hD hpd A,
    fun Cc _ γ => ⟨γ, (Category.id_comp γ).symm, fun β hβ => by
      have hg : γ = β := by simpa using hβ
      exact hg.symm⟩⟩

/-- ★**(i)(b)** 底の同型は pre-step の span に持ち上がる。

★`𝔽_Φ` の証人に `G` 成分 `1` を添えるだけ ——
★**結論の等式は `𝒟` の中の話なので `G` 成分に影響されない。** -/
theorem twPreStepSpan (A B : TwObj Φ G)
    (α : ((twPreFrobenioid (G := G) hD hpd).toElem.obj A).base ⟶
      ((twPreFrobenioid (G := G) hD hpd).toElem.obj B).base) (hα : IsIso α) :
    ∃ (X : TwObj Φ G) (φ : X ⟶ A) (ψ : X ⟶ B)
      (hφ : IsPreStep (twPreFrobenioid (G := G) hD hpd) φ),
      IsPreStep (twPreFrobenioid (G := G) hD hpd) ψ ∧
      α = @inv _ _ _ _ ((twPreFrobenioid (G := G) hD hpd).Base φ) hφ.2
            ≫ (twPreFrobenioid (G := G) hD hpd).Base ψ := by
  obtain ⟨X₀, φ₀, ψ₀, hφ₀, hψ₀, heq⟩ :=
    elemFrob_preStepSpan Φ hD hpd A.ofElem B.ofElem α hα
  exact ⟨⟨X₀⟩, ⟨φ₀, 1⟩, ⟨ψ₀, 1⟩, hφ₀, hψ₀, heq⟩

/-- ★★**pull-back 性は `𝔽_Φ` 成分だけで決まる**(次数 1 のとき)。

★**単射性**: `G` 成分は次数 1 なので消約される。
★**全射性**: 目標の `G` 成分から `f.unit⁻¹` を掛けて逆算する。 -/
theorem twIsPullBack_of {A B : TwObj Φ G} (f : A ⟶ B) (hlin : f.hom.deg = 1)
    (h : IsPullBack (elemPreFrobenioid Φ hD hpd) f.hom) :
    IsPullBack (twPreFrobenioid (G := G) hD hpd) f := by
  intro X
  obtain ⟨hinj, hsurj⟩ := h X.ofElem
  constructor
  · intro g₁ g₂ hg
    have hv := congrArg Subtype.val hg
    have h1 : g₁ ≫ f = g₂ ≫ f := congrArg Prod.fst hv
    have h2 : g₁.hom.base = g₂.hom.base := congrArg Prod.snd hv
    have hhom : g₁.hom = g₂.hom :=
      hinj (Subtype.ext (Prod.ext (congrArg TwHom.hom h1) h2))
    refine TwHom.ext hhom ?_
    have hu := congrArg TwHom.unit h1
    rw [twComp_unit, twComp_unit, hlin] at hu
    simpa using mul_right_cancel hu
  · intro y
    obtain ⟨f₀, hf₀⟩ := hsurj ⟨(y.1.1.hom, y.1.2), y.2⟩
    have hf₀' := congrArg Subtype.val hf₀
    have e1 : f₀ ≫ f.hom = y.1.1.hom := congrArg Prod.fst hf₀'
    have e2 : f₀.base = y.1.2 := congrArg Prod.snd hf₀'
    refine ⟨⟨f₀, y.1.1.unit * f.unit⁻¹⟩, Subtype.ext (Prod.ext ?_ e2)⟩
    refine TwHom.ext e1 ?_
    show (y.1.1.unit * f.unit⁻¹) ^ ((f.hom.deg : ℕ+) : ℕ) * f.unit = y.1.1.unit
    rw [hlin]
    simp

/-- ★★**pull-back 性の逆向きの橋** —— 捻れ積で pull-back なら `𝔽_Φ` でも pull-back。

★`G` 成分を `1` に取って `𝔽_Φ` の四角形を捻れ積に持ち上げるだけ。 -/
theorem twPullBack_hom {A B : TwObj Φ G} (α : A ⟶ B)
    (h : IsPullBack (twPreFrobenioid (G := G) hD hpd) α) :
    IsPullBack (elemPreFrobenioid Φ hD hpd) α.hom := by
  intro X₀
  obtain ⟨hinj, hsurj⟩ := h ⟨X₀⟩
  constructor
  · intro g₁ g₂ hg
    have hv := congrArg Subtype.val hg
    have e1 : g₁ ≫ α.hom = g₂ ≫ α.hom := congrArg Prod.fst hv
    have e2 : g₁.base = g₂.base := congrArg Prod.snd hv
    have := hinj (a₁ := (⟨g₁, 1⟩ : (⟨X₀⟩ : TwObj Φ G) ⟶ A))
      (a₂ := (⟨g₂, 1⟩ : (⟨X₀⟩ : TwObj Φ G) ⟶ A))
      (Subtype.ext (Prod.ext (TwHom.ext e1 rfl) e2))
    exact congrArg TwHom.hom this
  · intro y
    obtain ⟨g, hg⟩ := hsurj ⟨((⟨y.1.1, 1⟩ : (⟨X₀⟩ : TwObj Φ G) ⟶ B), y.1.2), y.2⟩
    have hg' := congrArg Subtype.val hg
    have hp1 : g ≫ α = (⟨y.1.1, 1⟩ : (⟨X₀⟩ : TwObj Φ G) ⟶ B) := congrArg Prod.fst hg'
    have hp2 : g.hom.base = y.1.2 := congrArg Prod.snd hg'
    exact ⟨g.hom, Subtype.ext (Prod.ext (congrArg TwHom.hom hp1) hp2)⟩

/-- ★★**pull-back 性の判定は `𝔽_Φ` と完全に一致する**。 -/
theorem twIsPullBack_iff {A B : TwObj Φ G} (f : A ⟶ B) :
    IsPullBack (twPreFrobenioid (G := G) hD hpd) f ↔
      ElemFrobCat.Hom.deg f.hom = 1 ∧ toChar (ElemFrobCat.Hom.div f.hom) = 0 := by
  constructor
  · intro h
    exact (elemFrob_isPullBack_iff Φ hD hpd f.hom).mp (twPullBack_hom hD hpd f h)
  · rintro ⟨hd, hz⟩
    exact twIsPullBack_of hD hpd f hd ((elemFrob_isPullBack_iff Φ hD hpd f.hom).mpr ⟨hd, hz⟩)

/-- ★**(iv)(b)** pull-back は LB-invertible かつ linear。

★co-angular は無条件、isometric と linear は `𝔽_Φ` 側から。 -/
theorem twPullBackLB {A B : TwObj Φ G} (α : A ⟶ B)
    (h : IsPullBack (twPreFrobenioid (G := G) hD hpd) α) :
    IsLBInvertible (twPreFrobenioid (G := G) hD hpd) α ∧
      IsLinear (twPreFrobenioid (G := G) hD hpd) α := by
  obtain ⟨hlb, hlin⟩ := elemFrob_pullBackLB Φ hD hpd α.hom (twPullBack_hom hD hpd α h)
  exact ⟨⟨twCoAngular hD hpd α, hlb.2⟩, hlin⟩

/-- ★**(iv)(a)** 任意射の 3 分解。

★`𝔽_Φ` の分解に `G` 成分を **Frobenius 因子へ集める** ——
残り 2 因子は次数 1 なので `G` 成分が素通りする。 -/
theorem twArbFactor {A B : TwObj Φ G} (φ : A ⟶ B) :
    ∃ (X Y : TwObj Φ G) (γ : A ⟶ X) (β : X ⟶ Y) (α : Y ⟶ B),
      φ = γ ≫ β ≫ α ∧ IsFrobeniusType (twPreFrobenioid (G := G) hD hpd) γ ∧
        IsPreStep (twPreFrobenioid (G := G) hD hpd) β ∧
        IsPullBack (twPreFrobenioid (G := G) hD hpd) α := by
  obtain ⟨X₀, Y₀, γ₀, β₀, α₀, hfac, hγ₀, hβ₀, hα₀⟩ := elemFrob_arbFactor Φ hD hpd φ.hom
  obtain ⟨-, hα₀lin⟩ := elemFrob_pullBackLB Φ hD hpd α₀ hα₀
  refine ⟨⟨X₀⟩, ⟨Y₀⟩, ⟨γ₀, φ.unit⟩, ⟨β₀, 1⟩, ⟨α₀, 1⟩, ?_,
    ⟨⟨twCoAngular hD hpd _, hγ₀.1.2⟩, hγ₀.2⟩, hβ₀,
    twIsPullBack_of hD hpd _ hα₀lin hα₀⟩
  refine TwHom.ext hfac ?_
  show φ.unit = φ.unit ^ ((ElemFrobCat.Hom.deg (β₀ ≫ α₀) : ℕ+) : ℕ) * ((1 : G) ^ _ * 1)
  rw [show ElemFrobCat.Hom.deg (β₀ ≫ α₀) = 1 by
    rw [ElemFrobCat.degFr_comp, show ElemFrobCat.Hom.deg α₀ = 1 from hα₀lin,
      show ElemFrobCat.Hom.deg β₀ = 1 from hβ₀.1, one_mul]]
  simp

/-- ★**(vi)** 単元を除く忠実性。

★`𝔽_Φ` は忠実なので `hom` は `α₀` の差だけ。
★**`G` 成分の差 `ψ.unit⁻¹ * φ.unit` を足せば単元になる。** -/
theorem twFaithfulUpToUnits {A B : TwObj Φ G} (φ ψ : A ⟶ B)
    (hb : BaseEquivalent (twPreFrobenioid (G := G) hD hpd) φ ψ)
    (hm : MetricallyEquivalent (twPreFrobenioid (G := G) hD hpd) φ ψ)
    (hφ : IsPreStep (twPreFrobenioid (G := G) hD hpd) φ)
    (hψ : IsPreStep (twPreFrobenioid (G := G) hD hpd) ψ) :
    ∃ α : End B, α ∈ OTimes (twPreFrobenioid (G := G) hD hpd) B ∧
      φ = ψ ≫ (α : B ⟶ B) := by
  obtain ⟨α₀, hα₀m, hα₀⟩ :=
    elemFrob_faithfulUpToUnits Φ hD hpd φ.hom ψ.hom hb hm hφ hψ
  haveI : IsIso (α₀ : B.ofElem ⟶ B.ofElem) :=
    (CategoryTheory.isUnit_iff_isIso (α₀ : End B.ofElem)).mp hα₀m.2
  refine ⟨⟨α₀, ψ.unit⁻¹ * φ.unit⟩, ⟨⟨hα₀m.1.1, hα₀m.1.2⟩, ?_⟩, ?_⟩
  · exact (CategoryTheory.isUnit_iff_isIso _).mpr (twIsIso_of _ ‹_›)
  · refine TwHom.ext hα₀ ?_
    show φ.unit = ψ.unit ^ ((ElemFrobCat.Hom.deg (α₀ : B.ofElem ⟶ B.ofElem) : ℕ+) : ℕ)
      * (ψ.unit⁻¹ * φ.unit)
    rw [show ElemFrobCat.Hom.deg (α₀ : B.ofElem ⟶ B.ofElem) = 1 from hα₀m.1.2]
    simp

/-- ★**(v)(b)** pre-step は「co-angular ∘ isometric」に分解する。

★`G` 成分は**第1因子へ集める**(第2因子は次数 1 なので素通り)。 -/
theorem twPreStepFactor {A B : TwObj Φ G} (φ : A ⟶ B)
    (hφ : IsPreStep (twPreFrobenioid (G := G) hD hpd) φ) :
    ∃ (X : TwObj Φ G) (β : A ⟶ X) (α : X ⟶ B),
      φ = β ≫ α ∧ IsCoAngular (twPreFrobenioid (G := G) hD hpd) β ∧
        IsPreStep (twPreFrobenioid (G := G) hD hpd) β ∧
        IsIsometric (twPreFrobenioid (G := G) hD hpd) α ∧
        IsPreStep (twPreFrobenioid (G := G) hD hpd) α := by
  obtain ⟨X₀, β₀, α₀, hfac, -, hβ₀s, hα₀i, hα₀s⟩ :=
    elemFrob_preStepFactor Φ hD hpd φ.hom hφ
  refine ⟨⟨X₀⟩, ⟨β₀, φ.unit⟩, ⟨α₀, 1⟩, ?_, twCoAngular hD hpd _, hβ₀s, hα₀i, hα₀s⟩
  refine TwHom.ext hfac ?_
  show φ.unit = φ.unit ^ ((ElemFrobCat.Hom.deg α₀ : ℕ+) : ℕ) * 1
  rw [show ElemFrobCat.Hom.deg α₀ = 1 from hα₀s.1]
  simp

/-- ★**(v)(c)** pre-step は「isometric ∘ co-angular」にも分解する。 -/
theorem twPreStepFactor' {A B : TwObj Φ G} (φ : A ⟶ B)
    (hφ : IsPreStep (twPreFrobenioid (G := G) hD hpd) φ) :
    ∃ (X : TwObj Φ G) (β : A ⟶ X) (α : X ⟶ B),
      φ = β ≫ α ∧ IsIsometric (twPreFrobenioid (G := G) hD hpd) β ∧
        IsPreStep (twPreFrobenioid (G := G) hD hpd) β ∧
        IsCoAngular (twPreFrobenioid (G := G) hD hpd) α ∧
        IsPreStep (twPreFrobenioid (G := G) hD hpd) α := by
  obtain ⟨X₀, β₀, α₀, hfac, hβ₀i, hβ₀s, -, hα₀s⟩ :=
    elemFrob_preStepFactor' Φ hD hpd φ.hom hφ
  refine ⟨⟨X₀⟩, ⟨β₀, φ.unit⟩, ⟨α₀, 1⟩, ?_, hβ₀i, hβ₀s, twCoAngular hD hpd _, hα₀s⟩
  refine TwHom.ext hfac ?_
  show φ.unit = φ.unit ^ ((ElemFrobCat.Hom.deg α₀ : ℕ+) : ℕ) * 1
  rw [show ElemFrobCat.Hom.deg α₀ = 1 from hα₀s.1]
  simp

/-! ### ★`𝒪^▷` の対応

★★**`OTri` の元は次数 1 なので、`G` 成分の指数がすべて消える。**
残るのは `G` の可換性だけで、対応は「`G` 成分をそのまま写す」になる。 -/

/-- ★`OTri` の帰属は `𝔽_Φ` 成分だけで決まる(定義の展開)。 -/
theorem tw_mem_otri_iff {A : TwObj Φ G} (α : End A) :
    α ∈ OTri (twPreFrobenioid (G := G) hD hpd) A ↔
      (α.hom : End A.ofElem) ∈ OTri (elemPreFrobenioid Φ hD hpd) A.ofElem := Iff.rfl

/-- ★**(iii)(c)** 順方向。★`G` 成分は `α.unit` をそのまま渡す。 -/
theorem twOtriFwd {A B : TwObj Φ G} (φ : A ⟶ B)
    (_ : IsCoAngular (twPreFrobenioid (G := G) hD hpd) φ)
    (hst : IsPreStep (twPreFrobenioid (G := G) hD hpd) φ)
    (α : End A) (hα : α ∈ OTri (twPreFrobenioid (G := G) hD hpd) A) :
    ∃! β : End B, β ∈ OTri (twPreFrobenioid (G := G) hD hpd) B ∧
      (φ ≫ β : A ⟶ B) = (α : A ⟶ A) ≫ φ := by
  have hφd : ElemFrobCat.Hom.deg φ.hom = 1 := hst.1
  have hαd : ElemFrobCat.Hom.deg α.hom = 1 := ((tw_mem_otri_iff hD hpd α).mp hα).2
  obtain ⟨β₀, ⟨hβ₀m, hβ₀e⟩, huniq⟩ :=
    elemFrob_otriFwd Φ hD hpd φ.hom hst α.hom ((tw_mem_otri_iff hD hpd α).mp hα)
  have hβ₀d : ElemFrobCat.Hom.deg β₀ = 1 := hβ₀m.2
  refine ⟨⟨β₀, α.unit⟩, ⟨hβ₀m, ?_⟩, ?_⟩
  · refine TwHom.ext hβ₀e ?_
    show φ.unit ^ ((ElemFrobCat.Hom.deg β₀ : ℕ+) : ℕ) * α.unit
      = α.unit ^ ((ElemFrobCat.Hom.deg φ.hom : ℕ+) : ℕ) * φ.unit
    rw [hβ₀d, hφd]
    simpa using mul_comm φ.unit α.unit
  · rintro β ⟨hβ, hβe⟩
    have hβd : ElemFrobCat.Hom.deg β.hom = 1 := ((tw_mem_otri_iff hD hpd β).mp hβ).2
    have hhom : β.hom = β₀ :=
      huniq β.hom ⟨(tw_mem_otri_iff hD hpd β).mp hβ, congrArg TwHom.hom hβe⟩
    refine TwHom.ext hhom ?_
    have hu := congrArg TwHom.unit hβe
    rw [twComp_unit, twComp_unit, hβd, hφd] at hu
    simp only [PNat.one_coe, pow_one] at hu
    exact mul_left_cancel (a := φ.unit) (by rw [hu, mul_comm])

/-- ★**(iii)(c)** 逆方向。★`G` 成分は `β.unit` をそのまま戻す。 -/
theorem twOtriBwd {A B : TwObj Φ G} (φ : A ⟶ B)
    (_ : IsCoAngular (twPreFrobenioid (G := G) hD hpd) φ)
    (hst : IsPreStep (twPreFrobenioid (G := G) hD hpd) φ)
    (β : End B) (hβ : β ∈ OTri (twPreFrobenioid (G := G) hD hpd) B) :
    ∃! α : End A, α ∈ OTri (twPreFrobenioid (G := G) hD hpd) A ∧
      (φ ≫ β : A ⟶ B) = (α : A ⟶ A) ≫ φ := by
  have hφd : ElemFrobCat.Hom.deg φ.hom = 1 := hst.1
  have hβd : ElemFrobCat.Hom.deg β.hom = 1 := ((tw_mem_otri_iff hD hpd β).mp hβ).2
  obtain ⟨α₀, ⟨hα₀m, hα₀e⟩, huniq⟩ :=
    elemFrob_otriBwd Φ hD hpd φ.hom hst β.hom ((tw_mem_otri_iff hD hpd β).mp hβ)
  have hα₀d : ElemFrobCat.Hom.deg α₀ = 1 := hα₀m.2
  refine ⟨⟨α₀, β.unit⟩, ⟨hα₀m, ?_⟩, ?_⟩
  · refine TwHom.ext hα₀e ?_
    show φ.unit ^ ((ElemFrobCat.Hom.deg β.hom : ℕ+) : ℕ) * β.unit
      = β.unit ^ ((ElemFrobCat.Hom.deg φ.hom : ℕ+) : ℕ) * φ.unit
    rw [hβd, hφd]
    simpa using mul_comm φ.unit β.unit
  · rintro α ⟨hα, hαe⟩
    have hαd : ElemFrobCat.Hom.deg α.hom = 1 := ((tw_mem_otri_iff hD hpd α).mp hα).2
    have hhom : α.hom = α₀ :=
      huniq α.hom ⟨(tw_mem_otri_iff hD hpd α).mp hα, congrArg TwHom.hom hαe⟩
    refine TwHom.ext hhom ?_
    have hu := congrArg TwHom.unit hαe
    rw [twComp_unit, twComp_unit, hβd, hφd] at hu
    simp only [PNat.one_coe, pow_one] at hu
    exact mul_right_cancel (b := φ.unit) (by rw [← hu, mul_comm])

/-- ★**(iii)(c)** 対応は `Base φ` にしか依らない。

★`G` 成分の等式は `φ.unit * β.unit = α.unit * φ.unit`、すなわち
★**`β.unit = α.unit`** であり、`φ` を `φ'` に取り替えても同じ形になる。 -/
theorem twOtriBase {A B : TwObj Φ G} (φ φ' : A ⟶ B)
    (_ : IsCoAngular (twPreFrobenioid (G := G) hD hpd) φ)
    (hst : IsPreStep (twPreFrobenioid (G := G) hD hpd) φ)
    (_ : IsCoAngular (twPreFrobenioid (G := G) hD hpd) φ')
    (hst' : IsPreStep (twPreFrobenioid (G := G) hD hpd) φ')
    (hbase : (twPreFrobenioid (G := G) hD hpd).Base φ
      = (twPreFrobenioid (G := G) hD hpd).Base φ')
    (α : End A) (hα : α ∈ OTri (twPreFrobenioid (G := G) hD hpd) A)
    (β : End B) (hβ : β ∈ OTri (twPreFrobenioid (G := G) hD hpd) B)
    (h : (φ ≫ β : A ⟶ B) = (α : A ⟶ A) ≫ φ) :
    (φ' ≫ β : A ⟶ B) = (α : A ⟶ A) ≫ φ' := by
  have hφd : ElemFrobCat.Hom.deg φ.hom = 1 := hst.1
  have hφd' : ElemFrobCat.Hom.deg φ'.hom = 1 := hst'.1
  have hαd : ElemFrobCat.Hom.deg α.hom = 1 := ((tw_mem_otri_iff hD hpd α).mp hα).2
  have hβd : ElemFrobCat.Hom.deg β.hom = 1 := ((tw_mem_otri_iff hD hpd β).mp hβ).2
  refine TwHom.ext (elemFrob_otriBase Φ hD hpd φ.hom φ'.hom hst hst' hbase
    α.hom ((tw_mem_otri_iff hD hpd α).mp hα)
    β.hom ((tw_mem_otri_iff hD hpd β).mp hβ) (congrArg TwHom.hom h)) ?_
  have hu := congrArg TwHom.unit h
  rw [twComp_unit, twComp_unit, hβd, hφd] at hu
  show φ'.unit ^ ((ElemFrobCat.Hom.deg β.hom : ℕ+) : ℕ) * β.unit
    = α.unit ^ ((ElemFrobCat.Hom.deg φ'.hom : ℕ+) : ℕ) * φ'.unit
  rw [hβd, hφd']
  simp only [PNat.one_coe, pow_one] at hu ⊢
  rw [show β.unit = α.unit from mul_left_cancel (a := φ.unit) (by rw [hu, mul_comm])]
  exact mul_comm _ _

/-! ### ★分解の一意性

★★**同型 `γ` の `G` 成分は「両辺の `G` 成分の差」で一意に決まる。**
2 つの条件(`α' = γ.inv ≫ α` と `β' = β ≫ γ.hom`)が同じ `g` を要求するが、
★**それが一致するのは仮定 `β ≫ α = β' ≫ α'` の `G` 成分そのものである。** -/

/-- ★`𝔽_Φ` の同型に `G` 成分を足して捻れ積の同型にする。 -/
def twIsoOf {X X' : ElemFrobCat Φ} (e : X ≅ X') (g : G) :
    (⟨X⟩ : TwObj Φ G) ≅ ⟨X'⟩ where
  hom := ⟨e.hom, g⟩
  inv := ⟨e.inv, g⁻¹⟩
  hom_inv_id := by
    have hdh : ElemFrobCat.Hom.deg e.hom = 1 :=
      ((ElemFrobCat.isIso_iff e.hom).mp inferInstance).2.2
    have hd : ElemFrobCat.Hom.deg e.inv = 1 := by
      have h1 := congrArg ElemFrobCat.Hom.deg e.hom_inv_id
      rw [ElemFrobCat.degFr_comp, hdh, mul_one] at h1
      simpa using h1
    refine TwHom.ext e.hom_inv_id ?_
    show g ^ ((ElemFrobCat.Hom.deg e.inv : ℕ+) : ℕ) * g⁻¹ = (1 : G)
    rw [hd]; simp
  inv_hom_id := by
    have hd : ElemFrobCat.Hom.deg e.hom = 1 :=
      ((ElemFrobCat.isIso_iff e.hom).mp inferInstance).2.2
    refine TwHom.ext e.inv_hom_id ?_
    show g⁻¹ ^ ((ElemFrobCat.Hom.deg e.hom : ℕ+) : ℕ) * g = (1 : G)
    rw [hd]; simp

@[simp] theorem twIsoOf_hom_unit {X X' : ElemFrobCat Φ} (e : X ≅ X') (g : G) :
    (twIsoOf (G := G) e g).hom.unit = g := rfl

@[simp] theorem twIsoOf_inv_unit {X X' : ElemFrobCat Φ} (e : X ≅ X') (g : G) :
    (twIsoOf (G := G) e g).inv.unit = g⁻¹ := rfl

/-- ★**(v)(b)** の一意性。 -/
theorem twPreStepFactorUniq {A B : TwObj Φ G} (X X' : TwObj Φ G)
    (β : A ⟶ X) (α : X ⟶ B) (β' : A ⟶ X') (α' : X' ⟶ B) (heq : β ≫ α = β' ≫ α')
    (_ : IsCoAngular (twPreFrobenioid (G := G) hD hpd) β)
    (hβs : IsPreStep (twPreFrobenioid (G := G) hD hpd) β)
    (hαi : IsIsometric (twPreFrobenioid (G := G) hD hpd) α)
    (hαs : IsPreStep (twPreFrobenioid (G := G) hD hpd) α)
    (_ : IsCoAngular (twPreFrobenioid (G := G) hD hpd) β')
    (hβs' : IsPreStep (twPreFrobenioid (G := G) hD hpd) β')
    (hαi' : IsIsometric (twPreFrobenioid (G := G) hD hpd) α')
    (hαs' : IsPreStep (twPreFrobenioid (G := G) hD hpd) α') :
    ∃ γ : X ≅ X', α' = γ.inv ≫ α ∧ β' = β ≫ γ.hom := by
  obtain ⟨γ₀, hα₀, hβ₀⟩ := elemFrob_preStepFactorUniq Φ hD hpd X.ofElem X'.ofElem
    β.hom α.hom β'.hom α'.hom (congrArg TwHom.hom heq)
    hβs hαi hαs.1 hβs' hαi' hαs'.1
  have hu := congrArg TwHom.unit heq
  rw [twComp_unit, twComp_unit, show ElemFrobCat.Hom.deg α.hom = 1 from hαs.1,
    show ElemFrobCat.Hom.deg α'.hom = 1 from hαs'.1] at hu
  simp only [PNat.one_coe, pow_one] at hu
  refine ⟨twIsoOf (G := G) γ₀ (β.unit⁻¹ * β'.unit), ?_, ?_⟩
  · refine TwHom.ext hα₀ ?_
    show α'.unit = (β.unit⁻¹ * β'.unit)⁻¹ ^ ((ElemFrobCat.Hom.deg α.hom : ℕ+) : ℕ) * α.unit
    rw [show ElemFrobCat.Hom.deg α.hom = 1 from hαs.1]
    simp only [PNat.one_coe, pow_one, mul_inv_rev, inv_inv]
    rw [mul_assoc, hu]
    simp
  · refine TwHom.ext hβ₀ ?_
    show β'.unit = β.unit ^ ((ElemFrobCat.Hom.deg γ₀.hom : ℕ+) : ℕ) * (β.unit⁻¹ * β'.unit)
    rw [show ElemFrobCat.Hom.deg γ₀.hom = 1 from
      ((ElemFrobCat.isIso_iff γ₀.hom).mp inferInstance).2.2]
    simp

/-- ★**(v)(c)** の一意性。★同じ計算(`α` 側も `β` 側も次数 1)。 -/
theorem twPreStepFactorUniq' {A B : TwObj Φ G} (X X' : TwObj Φ G)
    (β : A ⟶ X) (α : X ⟶ B) (β' : A ⟶ X') (α' : X' ⟶ B) (heq : β ≫ α = β' ≫ α')
    (hβi : IsIsometric (twPreFrobenioid (G := G) hD hpd) β)
    (hβs : IsPreStep (twPreFrobenioid (G := G) hD hpd) β)
    (_ : IsCoAngular (twPreFrobenioid (G := G) hD hpd) α)
    (hαs : IsPreStep (twPreFrobenioid (G := G) hD hpd) α)
    (hβi' : IsIsometric (twPreFrobenioid (G := G) hD hpd) β')
    (hβs' : IsPreStep (twPreFrobenioid (G := G) hD hpd) β')
    (_ : IsCoAngular (twPreFrobenioid (G := G) hD hpd) α')
    (hαs' : IsPreStep (twPreFrobenioid (G := G) hD hpd) α') :
    ∃ γ : X ≅ X', α' = γ.inv ≫ α ∧ β' = β ≫ γ.hom := by
  obtain ⟨γ₀, hα₀, hβ₀⟩ := elemFrob_preStepFactorUniq' Φ hD hpd X.ofElem X'.ofElem
    β.hom α.hom β'.hom α'.hom (congrArg TwHom.hom heq)
    hβi hβs hαs hβi' hβs' hαs'
  have hu := congrArg TwHom.unit heq
  rw [twComp_unit, twComp_unit, show ElemFrobCat.Hom.deg α.hom = 1 from hαs.1,
    show ElemFrobCat.Hom.deg α'.hom = 1 from hαs'.1] at hu
  simp only [PNat.one_coe, pow_one] at hu
  refine ⟨twIsoOf (G := G) γ₀ (β.unit⁻¹ * β'.unit), ?_, ?_⟩
  · refine TwHom.ext hα₀ ?_
    show α'.unit = (β.unit⁻¹ * β'.unit)⁻¹ ^ ((ElemFrobCat.Hom.deg α.hom : ℕ+) : ℕ) * α.unit
    rw [show ElemFrobCat.Hom.deg α.hom = 1 from hαs.1]
    simp only [PNat.one_coe, pow_one, mul_inv_rev, inv_inv]
    rw [mul_assoc, hu]
    simp
  · refine TwHom.ext hβ₀ ?_
    show β'.unit = β.unit ^ ((ElemFrobCat.Hom.deg γ₀.hom : ℕ+) : ℕ) * (β.unit⁻¹ * β'.unit)
    rw [show ElemFrobCat.Hom.deg γ₀.hom = 1 from
      ((ElemFrobCat.isIso_iff γ₀.hom).mp inferInstance).2.2]
    simp

/-- ★**(iv)(a)** 3 分解の一意性。

★★**`G` 成分は 2 つの同型が別々に決める** ——
`ε` は Frobenius 因子の差 `γ.unit⁻¹ * γ'.unit`、
`δ` は pull-back 因子の差 `α.unit * α'.unit⁻¹`。
★**中央の条件が両者と整合するのは、まさに仮定の `G` 成分
`γ.unit * β.unit * α.unit = γ'.unit * β'.unit * α'.unit` である。** -/
theorem twArbFactorUniq {A B : TwObj Φ G} (X Y X' Y' : TwObj Φ G)
    (γ : A ⟶ X) (β : X ⟶ Y) (α : Y ⟶ B) (γ' : A ⟶ X') (β' : X' ⟶ Y') (α' : Y' ⟶ B)
    (heq : (γ ≫ β ≫ α : A ⟶ B) = γ' ≫ β' ≫ α')
    (hγ : IsFrobeniusType (twPreFrobenioid (G := G) hD hpd) γ)
    (hβ : IsPreStep (twPreFrobenioid (G := G) hD hpd) β)
    (hα : IsPullBack (twPreFrobenioid (G := G) hD hpd) α)
    (hγ' : IsFrobeniusType (twPreFrobenioid (G := G) hD hpd) γ')
    (hβ' : IsPreStep (twPreFrobenioid (G := G) hD hpd) β')
    (hα' : IsPullBack (twPreFrobenioid (G := G) hD hpd) α') :
    ∃ (δ : Y ≅ Y') (ε : X ≅ X'),
      α' = δ.inv ≫ α ∧ β' = ε.inv ≫ β ≫ δ.hom ∧ γ' = γ ≫ ε.hom := by
  have hα0 : IsPullBack (elemPreFrobenioid Φ hD hpd) α.hom := twPullBack_hom hD hpd α hα
  have hα0' : IsPullBack (elemPreFrobenioid Φ hD hpd) α'.hom := twPullBack_hom hD hpd α' hα'
  have hαl : ElemFrobCat.Hom.deg α.hom = 1 := (elemFrob_pullBackLB Φ hD hpd α.hom hα0).2
  have hαl' : ElemFrobCat.Hom.deg α'.hom = 1 := (elemFrob_pullBackLB Φ hD hpd α'.hom hα0').2
  have hβl : ElemFrobCat.Hom.deg β.hom = 1 := hβ.1
  have hβl' : ElemFrobCat.Hom.deg β'.hom = 1 := hβ'.1
  obtain ⟨δ₀, ε₀, hδα, hδβ, hδγ⟩ := elemFrob_arbFactorUniq Φ hD hpd
    X.ofElem Y.ofElem X'.ofElem Y'.ofElem γ.hom β.hom α.hom γ'.hom β'.hom α'.hom
    (congrArg TwHom.hom heq) (twFrobType_hom hD hpd γ hγ) hβ hα0
    (twFrobType_hom hD hpd γ' hγ') hβ' hα0'
  have hδd : ElemFrobCat.Hom.deg δ₀.hom = 1 :=
    ((ElemFrobCat.isIso_iff δ₀.hom).mp inferInstance).2.2
  have hεd : ElemFrobCat.Hom.deg ε₀.hom = 1 :=
    ((ElemFrobCat.isIso_iff ε₀.hom).mp inferInstance).2.2
  -- ★仮定の `G` 成分(すべての次数が 1 なので指数が消える)
  have hu : γ.unit * (β.unit * α.unit) = γ'.unit * (β'.unit * α'.unit) := by
    have h := congrArg TwHom.unit heq
    rw [twComp_unit, twComp_unit, twComp_unit, twComp_unit, twComp_hom, twComp_hom,
      ElemFrobCat.degFr_comp, ElemFrobCat.degFr_comp,
      hαl, hαl', hβl, hβl'] at h
    simpa using h
  refine ⟨twIsoOf (G := G) δ₀ (α.unit * α'.unit⁻¹),
    twIsoOf (G := G) ε₀ (γ.unit⁻¹ * γ'.unit), ?_, ?_, ?_⟩
  · refine TwHom.ext hδα ?_
    show α'.unit = (α.unit * α'.unit⁻¹)⁻¹ ^ ((ElemFrobCat.Hom.deg α.hom : ℕ+) : ℕ) * α.unit
    rw [hαl]
    simp [mul_comm, mul_assoc, mul_left_comm]
  · refine TwHom.ext hδβ ?_
    show β'.unit = (γ.unit⁻¹ * γ'.unit)⁻¹
        ^ ((ElemFrobCat.Hom.deg (β.hom ≫ δ₀.hom) : ℕ+) : ℕ)
        * (β.unit ^ ((ElemFrobCat.Hom.deg δ₀.hom : ℕ+) : ℕ) * (α.unit * α'.unit⁻¹))
    rw [ElemFrobCat.degFr_comp, hβl, hδd, mul_one]
    have hb' : β'.unit = γ'.unit⁻¹ * (γ.unit * (β.unit * α.unit)) * α'.unit⁻¹ := by
      rw [hu]; simp
    rw [hb']
    simp [mul_inv_rev, mul_comm, mul_assoc, mul_left_comm]
  · refine TwHom.ext hδγ ?_
    show γ'.unit = γ.unit ^ ((ElemFrobCat.Hom.deg ε₀.hom : ℕ+) : ℕ) * (γ.unit⁻¹ * γ'.unit)
    rw [hεd]
    simp

/-- ★**(i)(c)** `(𝒞^pl-bk)_A → 𝒟_{A_𝒟}` は圏同値。

★★**`G` 成分は三角形 `f ≫ W.hom = Z.hom` が一意に決める** ——
`W.hom` は pull-back で次数 1 なので `f.unit = Z.unit * W.unit⁻¹` である。
★**忠実性の `G` 成分はこの消約、充満性の `G` 成分は
`IsPullBack` の全射性が捻れ積のままで与える。** -/
theorem twPlBkEquiv (A : TwObj Φ G) :
    (plBkOverFunctor (twPreFrobenioid (G := G) hD hpd) A).IsEquivalence := by
  haveI hfaith : (plBkOverFunctor (twPreFrobenioid (G := G) hD hpd) A).Faithful := by
    constructor
    intro Z W f g hfg
    have hb : ElemFrobCat.Hom.base f.left.hom.hom = ElemFrobCat.Hom.base g.left.hom.hom :=
      congrArg CommaMorphism.left hfg
    obtain ⟨hWl, -⟩ := (twIsPullBack_iff hD hpd W.hom.hom).mp W.hom.property
    obtain ⟨hfl, -⟩ := (twIsPullBack_iff hD hpd f.left.hom).mp f.left.property
    obtain ⟨hgl, -⟩ := (twIsPullBack_iff hD hpd g.left.hom).mp g.left.property
    have hwf : (f.left.hom ≫ W.hom.hom) = Z.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Over.w f)
    have hwg : (g.left.hom ≫ W.hom.hom) = Z.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Over.w g)
    letI := isCancelAdd_of_isIntegralMonoid _ ((hpd Z.left.obj.ofElem.base).1)
    have hdiv : ElemFrobCat.Hom.div f.left.hom.hom = ElemFrobCat.Hom.div g.left.hom.hom := by
      have h1 := congrArg (fun x : _ ⟶ _ => ElemFrobCat.Hom.div (TwHom.hom x))
        (hwf.trans hwg.symm)
      simp only [twComp_hom] at h1
      rw [ElemFrobCat.div_comp, ElemFrobCat.div_comp, hb, hWl] at h1
      simp only [PNat.one_coe, one_smul] at h1
      exact add_left_cancel h1
    have hunit : f.left.hom.unit = g.left.hom.unit := by
      have hu := congrArg TwHom.unit (hwf.trans hwg.symm)
      rw [twComp_unit, twComp_unit, hWl] at hu
      simp only [PNat.one_coe, pow_one] at hu
      exact mul_right_cancel hu
    exact Over.OverMorphism.ext (InducedWideCategory.Hom.ext
      (TwHom.ext (ElemFrobCat.Hom.ext hb hdiv (hfl.trans hgl.symm)) hunit))
  haveI hfull : (plBkOverFunctor (twPreFrobenioid (G := G) hD hpd) A).Full := by
    constructor
    intro Z W h
    obtain ⟨hWl, hWi⟩ := (twIsPullBack_iff hD hpd W.hom.hom).mp W.hom.property
    obtain ⟨hZl, hZi⟩ := (twIsPullBack_iff hD hpd Z.hom.hom).mp Z.hom.property
    obtain ⟨S, hS⟩ := toChar_eq_zero_iff.mp hZi
    obtain ⟨f₀, hf₀⟩ := (W.hom.property Z.left.obj).2 ⟨(Z.hom.hom, h.left), (Over.w h).symm⟩
    have hp := Subtype.ext_iff.mp hf₀
    have h1 : (f₀ ≫ W.hom.hom) = Z.hom.hom := congrArg Prod.fst hp
    have h2 : ElemFrobCat.Hom.base f₀.hom = h.left := congrArg Prod.snd hp
    have h1' : (f₀.hom ≫ W.hom.hom.hom) = Z.hom.hom.hom := congrArg TwHom.hom h1
    have hdeg : ElemFrobCat.Hom.deg f₀.hom = 1 := by
      have hh := congrArg ElemFrobCat.Hom.deg h1'
      rw [ElemFrobCat.comp_deg, hWl, one_mul, hZl] at hh
      exact hh
    have hd : Φ.map (ElemFrobCat.Hom.base f₀.hom) (ElemFrobCat.Hom.div W.hom.hom.hom)
        + ElemFrobCat.Hom.div f₀.hom = ElemFrobCat.Hom.div Z.hom.hom.hom := by
      have hh := congrArg ElemFrobCat.Hom.div h1'
      rw [ElemFrobCat.div_comp, hWl] at hh
      simpa using hh
    have hdivu : IsAddUnit (ElemFrobCat.Hom.div f₀.hom) := by
      refine ⟨⟨ElemFrobCat.Hom.div f₀.hom,
        Φ.map (ElemFrobCat.Hom.base f₀.hom) (ElemFrobCat.Hom.div W.hom.hom.hom) + S.neg,
        ?_, ?_⟩, rfl⟩
      · rw [← add_assoc, add_comm (ElemFrobCat.Hom.div f₀.hom), hd, ← hS, S.val_neg]
      · rw [add_assoc, add_comm S.neg, ← add_assoc, hd, ← hS, S.val_neg]
    refine ⟨Over.homMk (⟨f₀, (twIsPullBack_iff hD hpd f₀).mpr
      ⟨hdeg, toChar_eq_zero_iff.mpr hdivu⟩⟩ : Z.left ⟶ W.left)
      (InducedWideCategory.Hom.ext h1), Over.OverMorphism.ext h2⟩
  haveI hess : (plBkOverFunctor (twPreFrobenioid (G := G) hD hpd) A).EssSurj := by
    constructor
    intro Y
    obtain ⟨q, hq⟩ : ∃ q : Y.left ⟶ A.ofElem.base, q = Y.hom := ⟨Y.hom, rfl⟩
    refine ⟨Over.mk (show (⟨(⟨⟨Y.left⟩⟩ : TwObj Φ G)⟩ :
        PlBk (twPreFrobenioid (G := G) hD hpd))
          ⟶ (⟨A⟩ : PlBk (twPreFrobenioid (G := G) hD hpd)) from
      ⟨(⟨(⟨q, 0, 1⟩ : (⟨Y.left⟩ : ElemFrobCat Φ) ⟶ A.ofElem), 1⟩ :
          (⟨⟨Y.left⟩⟩ : TwObj Φ G) ⟶ A),
        (twIsPullBack_iff hD hpd _).mpr ⟨rfl, map_zero _⟩⟩), ⟨?_⟩⟩
    refine Over.isoMk (Iso.refl _) ?_
    show 𝟙 Y.left ≫ Y.hom = q
    rw [Category.id_comp, hq]
  exact ⟨hfaith, hfull, hess⟩

/-! ## ★★★捻れ積は `Definition 1.3` の core 21 条をすべて満たす -/

/-- ★★★**`𝔽_Φ ⋉ G` は `FrobenioidCore` を満たす** —— 21 条すべて。 -/
theorem twFrobenioidCore : FrobenioidCore (twPreFrobenioid (G := G) hD hpd) where
  baseSurj := twBaseSurj hD hpd
  preStepSpan := twPreStepSpan hD hpd
  plBkEquiv := twPlBkEquiv hD hpd
  frobDegSurj := twFrobDegSurj hD hpd
  frobDegUniq := twFrobDegUniq hD hpd
  coAngularComp := twCoAngularComp hD hpd
  coAngularOfPreStep := fun α hca hst φ => twCoAngularOfPreStep hD hpd α hca hst φ
  otriFwd := fun φ hca hst α hα => twOtriFwd hD hpd φ hca hst α hα
  otriBwd := fun φ hca hst β hβ => twOtriBwd hD hpd φ hca hst β hβ
  otriBase := fun φ φ' hca hst hca' hst' hbase α hα β hβ h =>
    twOtriBase hD hpd φ φ' hca hst hca' hst' hbase α hα β hβ h
  arbFactor := twArbFactor hD hpd
  arbFactorUniq := fun X Y X' Y' γ β α γ' β' α' heq hγ hβ hα hγ' hβ' hα' =>
    twArbFactorUniq hD hpd X Y X' Y' γ β α γ' β' α' heq hγ hβ hα hγ' hβ' hα'
  pullBackLB := twPullBackLB hD hpd
  preStepMono := twPreStepMono hD hpd
  preStepFactor := twPreStepFactor hD hpd
  preStepFactorUniq := fun X X' β α β' α' heq hca hβ hαi hα hca' hβ' hαi' hα' =>
    twPreStepFactorUniq hD hpd X X' β α β' α' heq hca hβ hαi hα hca' hβ' hαi' hα'
  preStepFactor' := twPreStepFactor' hD hpd
  preStepFactorUniq' := fun X X' β α β' α' heq hβi hβ hca hα hβi' hβ' hca' hα' =>
    twPreStepFactorUniq' hD hpd X X' β α β' α' heq hβi hβ hca hα hβi' hβ' hca' hα'
  faithfulUpToUnits := fun φ ψ hb hm _ hφ _ hψ =>
    twFaithfulUpToUnits hD hpd φ ψ hb hm hφ hψ
  isotropicHullExists := twIsotropicHullExists hD hpd
  isotropicClosed := twIsotropicClosed hD hpd

/-! ## ★**(iii)(d)** —— 残る 2 本の圏同値

★★**目標は前順序(`OrderCat`)なので、示すことは「射が一意に決まる」ことだけ。**
捻れ積では三角形が `G` 成分まで一意に決めるので、★**前順序のまま**である。 -/

section CoaPre

variable [MorphismProperty.IsMultiplicative (coaPreProp (twPreFrobenioid (G := G) hD hpd))]

/-- ★**(iii)(d)** コスライス側。

★**忠実性は捻れ積が totally epimorphic であることから直接出る**
(`G` 成分の消約は `twTotallyEpimorphic` の中で済んでいる)。
★**充満性は `𝔽_Φ` 側の充満性に `G` 成分 `Z.unit⁻¹ * W.unit` を足すだけ。** -/
theorem twCoaPreUnderEquiv (A : TwObj Φ G) :
    (coaPreUnderFunctor (twPreFrobenioid (G := G) hD hpd) A).IsEquivalence := by
  haveI := coaPreProp_isMultiplicative (elemPreFrobenioid Φ hD hpd)
    (elemFrob_frobenioidCore Φ hD hpd).coAngularComp
  haveI := elemFrob_coaPreUnderEquiv Φ hD hpd A.ofElem
  haveI hfaith : (coaPreUnderFunctor (twPreFrobenioid (G := G) hD hpd) A).Faithful := by
    constructor
    intro Z W f g _
    have h1 : Z.hom.hom ≫ f.right.hom = W.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Under.w f)
    have h2 : Z.hom.hom ≫ g.right.hom = W.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Under.w g)
    haveI : Epi Z.hom.hom := twTotallyEpimorphic (G := G) hD (fun X => (hpd X).1) _ _ _
    exact Under.UnderMorphism.ext (InducedWideCategory.Hom.ext
      ((cancel_epi Z.hom.hom).mp (h1.trans h2.symm)))
  haveI hfull : (coaPreUnderFunctor (twPreFrobenioid (G := G) hD hpd) A).Full := by
    constructor
    intro Z W h
    -- ★`𝔽_Φ` 側の同じ図式(`G` 成分を落とす)
    obtain ⟨f₀, -⟩ := (coaPreUnderFunctor (elemPreFrobenioid Φ hD hpd) A.ofElem).map_surjective
      (X := Under.mk (show (⟨A.ofElem⟩ :
            WideSubcategory (coaPreProp (elemPreFrobenioid Φ hD hpd)))
          ⟶ (⟨Z.right.obj.ofElem⟩ :
            WideSubcategory (coaPreProp (elemPreFrobenioid Φ hD hpd))) from
        ⟨Z.hom.hom.hom, elemFrob_coAngular Φ hD hpd Z.hom.hom.hom,
          (Z.hom.property.2 : IsPreStep (elemPreFrobenioid Φ hD hpd) Z.hom.hom.hom)⟩))
      (Y := Under.mk (show (⟨A.ofElem⟩ :
            WideSubcategory (coaPreProp (elemPreFrobenioid Φ hD hpd)))
          ⟶ (⟨W.right.obj.ofElem⟩ :
            WideSubcategory (coaPreProp (elemPreFrobenioid Φ hD hpd))) from
        ⟨W.hom.hom.hom, elemFrob_coAngular Φ hD hpd W.hom.hom.hom,
          (W.hom.property.2 : IsPreStep (elemPreFrobenioid Φ hD hpd) W.hom.hom.hom)⟩))
      h
    have htri : Z.hom.hom.hom ≫ f₀.right.hom = W.hom.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Under.w f₀)
    haveI := Preorder.subsingleton_hom
      ((coaPreUnderFunctor (twPreFrobenioid (G := G) hD hpd) A).obj Z)
      ((coaPreUnderFunctor (twPreFrobenioid (G := G) hD hpd) A).obj W)
    refine ⟨Under.homMk (⟨(⟨f₀.right.hom, Z.hom.hom.unit⁻¹ * W.hom.hom.unit⟩ :
        Z.right.obj ⟶ W.right.obj),
      ⟨twCoAngular hD hpd _, f₀.right.property.2⟩⟩ : Z.right ⟶ W.right) ?_,
      Subsingleton.elim _ _⟩
    refine InducedWideCategory.Hom.ext ?_
    simp only [WideSubcategory.comp_def]
    refine TwHom.ext htri ?_
    show Z.hom.hom.unit ^ ((ElemFrobCat.Hom.deg f₀.right.hom : ℕ+) : ℕ)
      * (Z.hom.hom.unit⁻¹ * W.hom.hom.unit) = W.hom.hom.unit
    rw [show ElemFrobCat.Hom.deg f₀.right.hom = 1 from f₀.right.property.2.1]
    simp
  haveI hess : (coaPreUnderFunctor (twPreFrobenioid (G := G) hD hpd) A).EssSurj := by
    constructor
    intro c
    obtain ⟨a, ha⟩ : ∃ a : Φ.val A.ofElem.base, toChar a = c := by
      obtain ⟨y, hy⟩ := toChar_surjective _ c
      exact ⟨y, hy⟩
    refine ⟨Under.mk (show (⟨A⟩ : WideSubcategory (coaPreProp (twPreFrobenioid
        (G := G) hD hpd)))
        ⟶ (⟨A⟩ : WideSubcategory (coaPreProp (twPreFrobenioid (G := G) hD hpd))) from
      ⟨(⟨otriOf Φ A.ofElem a, 1⟩ : A ⟶ A), twCoAngular hD hpd _, rfl, ?_⟩), ⟨eqToIso ?_⟩⟩
    · show IsIso (𝟙 A.ofElem.base)
      infer_instance
    · exact congrArg toOrderCat ha
  exact ⟨hfaith, hfull, hess⟩

/-- ★**(iii)(d)** スライス側。

★**忠実性は `twPreStepMono`**(pre-step は次数 1 なので `G` 成分も消約できる)。 -/
theorem twCoaPreOverEquiv (A : TwObj Φ G) :
    (coaPreOverFunctor (twPreFrobenioid (G := G) hD hpd) A).IsEquivalence := by
  haveI := coaPreProp_isMultiplicative (elemPreFrobenioid Φ hD hpd)
    (elemFrob_frobenioidCore Φ hD hpd).coAngularComp
  haveI := elemFrob_coaPreOverEquiv Φ hD hpd A.ofElem
  haveI hfaith : (coaPreOverFunctor (twPreFrobenioid (G := G) hD hpd) A).Faithful := by
    constructor
    intro Z W f g _
    have h1 : f.left.hom ≫ W.hom.hom = Z.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Over.w f)
    have h2 : g.left.hom ≫ W.hom.hom = Z.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Over.w g)
    haveI : Mono W.hom.hom := twPreStepMono hD hpd _ W.hom.property.2
    exact Over.OverMorphism.ext (InducedWideCategory.Hom.ext
      ((cancel_mono W.hom.hom).mp (h1.trans h2.symm)))
  haveI hfull : (coaPreOverFunctor (twPreFrobenioid (G := G) hD hpd) A).Full := by
    constructor
    intro Z W h
    obtain ⟨f₀, -⟩ := (coaPreOverFunctor (elemPreFrobenioid Φ hD hpd) A.ofElem).map_surjective
      (X := Over.mk (show (⟨Z.left.obj.ofElem⟩ :
            WideSubcategory (coaPreProp (elemPreFrobenioid Φ hD hpd)))
          ⟶ (⟨A.ofElem⟩ :
            WideSubcategory (coaPreProp (elemPreFrobenioid Φ hD hpd))) from
        ⟨Z.hom.hom.hom, elemFrob_coAngular Φ hD hpd Z.hom.hom.hom,
          (Z.hom.property.2 : IsPreStep (elemPreFrobenioid Φ hD hpd) Z.hom.hom.hom)⟩))
      (Y := Over.mk (show (⟨W.left.obj.ofElem⟩ :
            WideSubcategory (coaPreProp (elemPreFrobenioid Φ hD hpd)))
          ⟶ (⟨A.ofElem⟩ :
            WideSubcategory (coaPreProp (elemPreFrobenioid Φ hD hpd))) from
        ⟨W.hom.hom.hom, elemFrob_coAngular Φ hD hpd W.hom.hom.hom,
          (W.hom.property.2 : IsPreStep (elemPreFrobenioid Φ hD hpd) W.hom.hom.hom)⟩))
      h
    have htri : f₀.left.hom ≫ W.hom.hom.hom = Z.hom.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Over.w f₀)
    haveI := Preorder.subsingleton_hom
      ((coaPreOverFunctor (twPreFrobenioid (G := G) hD hpd) A).obj W).unop
      ((coaPreOverFunctor (twPreFrobenioid (G := G) hD hpd) A).obj Z).unop
    refine ⟨Over.homMk (⟨(⟨f₀.left.hom, Z.hom.hom.unit * W.hom.hom.unit⁻¹⟩ :
        Z.left.obj ⟶ W.left.obj),
      ⟨twCoAngular hD hpd _, f₀.left.property.2⟩⟩ : Z.left ⟶ W.left) ?_,
      Quiver.Hom.unop_inj (Subsingleton.elim _ _)⟩
    refine InducedWideCategory.Hom.ext ?_
    simp only [WideSubcategory.comp_def]
    refine TwHom.ext htri ?_
    show (Z.hom.hom.unit * W.hom.hom.unit⁻¹)
      ^ ((ElemFrobCat.Hom.deg W.hom.hom.hom : ℕ+) : ℕ) * W.hom.hom.unit = Z.hom.hom.unit
    rw [show ElemFrobCat.Hom.deg W.hom.hom.hom = 1 from W.hom.property.2.1]
    simp
  haveI hess : (coaPreOverFunctor (twPreFrobenioid (G := G) hD hpd) A).EssSurj := by
    constructor
    intro c
    obtain ⟨a, ha⟩ : ∃ a : Φ.val A.ofElem.base, toChar a = c.unop := by
      obtain ⟨y, hy⟩ := toChar_surjective _ c.unop
      exact ⟨y, hy⟩
    refine ⟨Over.mk (show (⟨A⟩ : WideSubcategory (coaPreProp (twPreFrobenioid
        (G := G) hD hpd)))
        ⟶ (⟨A⟩ : WideSubcategory (coaPreProp (twPreFrobenioid (G := G) hD hpd))) from
      ⟨(⟨otriOf Φ A.ofElem a, 1⟩ : A ⟶ A), twCoAngular hD hpd _, rfl, ?_⟩), ⟨eqToIso ?_⟩⟩
    · show IsIso (𝟙 A.ofElem.base)
      infer_instance
    · refine Opposite.unop_injective ?_
      rw [← ha]
      exact elemFrob_map_inv_otriOf Φ hD hpd a _
  exact ⟨hfaith, hfull, hess⟩

end CoaPre

/-- ★★★**捻れ積 `𝔽_Φ ⋉ G` は Frobenioid である**(`Definition 1.3` の全条件)。 -/
theorem twIsFrobenioid : Frobenioid (twPreFrobenioid (G := G) hD hpd) := by
  haveI := coaPreProp_isMultiplicative (twPreFrobenioid (G := G) hD hpd)
    (twFrobenioidCore (G := G) hD hpd).coAngularComp
  exact ⟨twFrobenioidCore hD hpd, twCoaPreUnderEquiv hD hpd, twCoaPreOverEquiv hD hpd⟩

end Bridge

/-! ## ★★★具体的な反例 —— `𝔽_{ℕ on Vee} ⋉ (∏_n ℤ/n)`

★ここまでは `Φ`・`G` について一般だった。★**ここで両方を具体化して、
`Definition 1.3` をすべて満たす圏で「Frobenius 型射が mono でない」ことを出す。**

`𝒟 = Vee`(3点の半順序)、`Φ = ℕ`(定数)は witness で使っているものそのまま。

★★**`G` は `∏_{n} ℤ/n` を取る** —— ★**すべての次数 `d ≥ 2` について
`d`-捻れがある**ので、★★**次数 > 1 の射は 1 本も mono にならない。**
`ℤ/2` だけでは「次数 3 の射を使えばよい」という逃げ道が残るが、
★**この `G` はその逃げ道を塞ぐ。**

★原文 (FrdI p.42) の (a) の場合の議論はまさに
「次数をいくらでも大きくした prime-Frobenius 射 `ψ` を後置する」ものであり、
★**その `ψ` が FSMI(したがって mono)であることを要求している。** -/

/-- ★反例に使う `G` —— **すべての位数の捻れを持つ**可換群。 -/
abbrev CxG : Type := Multiplicative (∀ n : ℕ, ZMod n)

/-- ★★**任意の `d ≥ 2` について `d`-捻れ元がある**。 -/
theorem cxG_torsion (d : ℕ) (hd : 2 ≤ d) : ∃ u : CxG, u ≠ 1 ∧ u ^ d = 1 := by
  haveI : Fact (1 < d) := ⟨hd⟩
  refine ⟨Multiplicative.ofAdd (Pi.single d (1 : ZMod d)), ?_, ?_⟩
  · intro h
    have := congrFun (Multiplicative.ofAdd.injective h) d
    rw [Pi.single_eq_same] at this
    exact one_ne_zero this
  · refine Multiplicative.ofAdd.injective ?_
    show d • Pi.single d (1 : ZMod d) = 0
    funext n
    by_cases hn : n = d
    · subst hn
      simp only [Pi.smul_apply, Pi.single_eq_same, Pi.zero_apply, nsmul_eq_mul, mul_one]
      exact ZMod.natCast_self n
    · simp [Pi.single_eq_of_ne hn]

/-- ★反例の圏 `𝔽_{ℕ on Vee} ⋉ (∏_n ℤ/n)`。 -/
abbrev CxC : Type := TwObj wΦ CxG

/-- ★反例の pre-Frobenioid。 -/
def cxP : PreFrobenioid CxC wΦ.charOn :=
  twPreFrobenioid (G := CxG) isTotallyEpimorphic_vee (fun _ => isPreDivisorial_nat)

/-- ★★★**これは `Definition 1.3` の意味で Frobenioid である**。 -/
theorem cx_isFrobenioid : Frobenioid cxP :=
  twIsFrobenioid (G := CxG) isTotallyEpimorphic_vee (fun _ => isPreDivisorial_nat)

/-- ★★★**次数 > 1 の射は 1 本も mono でない**。

★★**これで「次数を大きくすれば mono になる」逃げ道が塞がる。** -/
theorem cx_not_mono_of_deg_gt_one {A B : CxC} (ζ : A ⟶ B) (hζ : ζ.unit = 1)
    (hd : 2 ≤ ((ζ.hom.deg : ℕ+) : ℕ)) : ¬ Mono ζ := by
  obtain ⟨u, hu, hup⟩ := cxG_torsion _ hd
  exact not_mono_twist ζ hζ u hu hup

/-- ★反例の対象。 -/
abbrev cxA : CxC := ⟨⟨Vee.top⟩⟩

/-- ★★★**次数 `n` の Frobenius 型自己射** —— `Definition 1.3, (ii)` が要求するもの。 -/
def cxZeta (n : ℕ+) : cxA ⟶ cxA := ⟨⟨𝟙 Vee.top, (0 : ℕ), n⟩, 1⟩

theorem cx_zeta_frobType (n : ℕ+) : IsFrobeniusType cxP (cxZeta n) := by
  refine ⟨⟨twCoAngular _ _ _, ?_⟩, ?_⟩
  · show toChar (0 : ℕ) = 0
    simp
  · show IsIso (𝟙 Vee.top)
    infer_instance

theorem cx_zeta_deg (n : ℕ+) : cxP.degFr (cxZeta n) = n := rfl

/-- ★★★**どの次数 > 1 でも mono にならない**。 -/
theorem cx_zeta_not_mono (n : ℕ+) (hn : 2 ≤ (n : ℕ)) : ¬ Mono (cxZeta n) :=
  cx_not_mono_of_deg_gt_one (cxZeta n) rfl hn

/-- ★★★**`Definition 1.3` は「Frobenius 型射は mono」を含意しない** ——
★**しかも次数をどれだけ大きくしても含意しない。**

★★**これが `Proposition 1.14, (iii)` の `⟸` を止めていたものの正体である。**
原文 (FrdI p.42) は (a) の場合に
「次数をいくらでも大きくした prime-Frobenius 射を後置する」と述べるが、
条件文はその射自身にも FSMI(したがって mono)を要求しており、
★`Definition 1.3` で mono を与える条項は `preStepMono` **1 つだけ**で、
pre-step にしか効かない。 -/
theorem cx_frobType_not_mono (n : ℕ+) (hn : 2 ≤ (n : ℕ)) :
    ∃ (A : CxC) (ζ : A ⟶ A),
      IsFrobeniusType cxP ζ ∧ cxP.degFr ζ = n ∧ ¬ Mono ζ :=
  ⟨cxA, cxZeta n, cx_zeta_frobType n, cx_zeta_deg n, cx_zeta_not_mono n hn⟩

/-! ### ★測定 —— ここまでで確定したこと(2026-08-16)

★**確定した**:
- 捻れ積は**圏である**(結合律・単位律を証明した)
- 射影関手 `twProj : 𝔽_Φ ⋉ G ⥤ 𝔽_Φ` がある
- ★★**次数 `d` の射(`G` 成分が `1`)は、`G` に `d`-捻れがあれば mono でない**
- ★**totally epimorphic である**(`G` 成分は群の消約)
- ★**connected である**(`𝔽_Φ` の証明をそのまま写せる)
- ★★**`PreFrobenioid` 構造を持つ**(`twPreFrobenioid`)
- ★★★**`Definition 1.3` の全条件を満たす**(`twFrobenioidCore` ＋ `twIsFrobenioid`)
- ★★★**具体化した `𝔽_{ℕ on Vee} ⋉ ℤ/2` で、次数 2 の Frobenius 型射が
  mono でない**(`cx_frobType_not_mono`)

★**したがって `Gap_1_14_iii` は ③(`sourceGap`)に上がる。**

★★**`sorry` は置かない。** 未確定のものは**書かない**。
-/

end ABC3.Check.FrdI
