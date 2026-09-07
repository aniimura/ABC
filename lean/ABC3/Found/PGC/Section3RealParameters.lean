import ABC3.Found.PGC.ArtinMap
import ABC3.Found.PGC.AbsClosureModules
import ABC3.Skeleton.PGC.Section3Defs
import Mathlib.NumberTheory.Cyclotomic.CyclotomicCharacter
import Mathlib.RepresentationTheory.Basic

/-!
# [pGC] §3 の自由パラメータを実物に固定する

`Skeleton/PGC/Section3.lean` の 2 定理は、原典に無い「項 `K` の自由な関数」を
仮説に取っている:

| 定理 | 自由なパラメータ | 型 |
|---|---|---|
| `cor_3_1` | `isHodgeTate` | `∀ K V …, Prop` |
| `cor_3_3` | `toGal` | `∀ K, {x : K.carrier // ‖x‖ = 1} → K.absGal` |

★どちらも `Check/PGC/FreeTermFunctionRefutation.lean`(2026-09-06)が
`sorry` 無しで**偽**だと示している(`not_cor_3_1_current_form` /
`not_cor_3_3_current_form`)。自由な項関数は同型な 2 つの項
(`selfField p` と `twistedField p`)に別の値を割り当ててよいからである。

本ファイルは**その 2 つのパラメータを実物に置き換える**。
`Check/PGC/Prop12ForallRD.lean` が `ResidueCardinality` について出した判断
(「量化を外して実物に固定せよ」)を §3 に対して実行したものにあたる。

## 1. `toGal` ← Artin 写像(§2–§6)

原文 (pGC p.6, Definition 3.2):
> We shall call the E[Γ_K]-module V uniformizing if the restriction of ρ_V to some open
> subgroup I of U_K (⊆ Γ^a_K^b) is the morphism I → E× induced by restricting some morphism
> of fields K → E to I ⊆ U_K ⊆ K.

★原典の `U_K ⊆ Γ_K^ab` は局所類体論の相互律を暗黙に使っている。
`Skeleton/PGC/Section3Defs.lean` はそこを「未構築の辞書 `toGal`」で代用していたが、
辞書はもう在庫にある —— `Found/PGC/ArtinMap.lean::exists_artinMap`(仮定ゼロ)。

### ★行き先の違い(自分で確かめた点)

`artinMap` の行き先は `Gal(K^ab/K) = (abelianClosure K ≃ₐ[K.carrier] abelianClosure K)`
であって `K.absGal` では**ない**。`IsUniformizing` は `ρ : K.absGal →* Eˣ` に
`toGal` を食わせるので、`Gal(K^ab/K)` の上では完結**しない**。
そこで `Γ_K ↠ Gal(K^ab/K)` の**全射性**(mathlib
`AlgEquiv.restrictNormalHom_surjective`、`K^ab/K` が正規であることは
`Found/PGC/AbelianClosure.lean::instNormalAbelianClosure`)を使って持ち上げる
——`artinToGal`。持ち上げは `Function.surjInv` による選択なので**一意ではない**が、
`IsUniformizing` が見るのは合成 `ρ ∘ toGal` だけなので、
`ρ` の側を `Gal(K^ab/K)` 経由の指標(`artinUnitChar`)に取れば選択に依らない
(`artinUnitChar_artinToGal`)。★これが「持ち上げが要るか」への答である:
**要る。ただし合成の値は選択に依らない。**

### 何が言えたか

* `isUniformizing_artin` —— 実物の `toGal` と実物の指標について
  Definition 3.2 が**真になる**(仮定ゼロ、`sorry` 無し)。
  開部分群は `I = U_K` そのもの、体準同型は `ι = id : K → K`。
* `artinToGal_injective` / `artinToGal_ne_const` —— 実物の `toGal` は単射であり、
  したがって**定数写像ではない**。★`Check/PGC/FreeTermFunctionRefutation.lean` の
  反例 `badToGal` は捻り側で定数 `1` を返すものだったから、
  固定すればあの反例は作れない。

## 2. `isHodgeTate` ← `CompKbar`(§7–§10)

原文 (pGC p.6, Corollary 3.1):
> Given a continuous Q_p[Γ_K]-vector space of finite Q_p-dimension, the issue of whether or
> not V is Hodge-Tate (as well as the invariants d_V(i)) can be determined entirely
> group-theoretically from the filtered group Γ_K.

原典の `d_V(i)` は `V(−i) ⊗_{Q_p} K̄^` の `Γ_K`-不変元の `K` 上の次元である。
Tate 捻り `V(−i)` を型として作る代わりに、**固有空間**として書いた:

    d_V(i) = dim_K { z ∈ ℂ_K ⊗_{Q_p} V | ∀ σ, (σ ⊗ σ) z = χ(σ)^i · z }

(`V(−i) ⊗ ℂ_K` の `Γ_K`-不変元 ⟺ `ℂ_K ⊗ V` の `χ^i`-固有元。型の同義語を作らずに済む。)
`ℂ_K` は `CompKbar K = closureCompletion K`(`Found/PGC/AbsClosureModules.lean`)。

### ★不足していた 1 つの instance(実測)

`Module ℚ_[p] (CompKbar K)` は**合成できなかった**:

```
error: failed to synthesize instance of type class
  Module ℚ_[p] (closureCompletion K)
```

原因は `UniformSpace.Completion.instModule`
(`Topology/Algebra/GroupCompletion.lean:204`)が `[UniformContinuousConstSMul R α]` を
要求していることで、`UniformContinuousConstSMul ℚ_[p] K.closure` は
★**`inferInstance` が失敗した(実測、`node tools/leanfile.mjs`)**。
一方 `Algebra ℚ_[p] K.closure`・`Module ℚ_[p] K.closure`・
`NormedAlgebra K.carrier K.closure`・`IsScalarTower ℚ_[p] K.carrier K.closure` は
どれも `inferInstance` が通る。そこで `closureNormedAlgebraQp`
(既存の `Algebra` を再利用して `norm_smul_le` だけ足す)を 1 つ立てたところ、
`UniformContinuousConstSMul` / `Module ℚ_[p] (CompKbar K)` /
`Algebra ℚ_[p] (CompKbar K)` / `IsScalarTower ℚ_[p] K.carrier (CompKbar K)` が
★**4 つとも `inferInstance` で通るようになった(実測)**
(繋いでいるのは `IsBoundedSMul.toUniformContinuousConstSMul`、
`Topology/MetricSpace/Algebra.lean:157` と見られる)。
★既存の `Algebra ℚ_[p] K.closure` を再利用しているので菱形は作っていない
——新しく `Algebra` を立てると `SMul` が
`UniformSpace.Completion.instSMul` と割れる(`lean-idioms.md` #263 に逐語で記録)。

### 何が言えたか / 言えていないか

* `hodgeTateWeightSpace` / `hodgeTateDim` / `IsHodgeTate` は**書ける**。
  自由パラメータは 1 つも残らない。
* `weightSpace_zero_ne_bot` —— 自明表現の重み 0 の空間は `0` でない(`1 ⊗ 1` が入る)。
* ★★**`d_V(i)` の値については何も言えていない。** `d_{triv}(0) = 1` すら
  `ℂ_K^{Γ_K} = K`(Ax–Sen–Tate)を要し、それは mathlib にも木にも無い。
  ★`Module.finrank` は無限次元で `0` を返すので、有限次元性の証明抜きでは
  `d_V(i) ≥ 1` すら出ない。

## 逸脱の記録(CLAUDE.md「逸脱」)

1. **`toGal` の持ち上げは選択を含む**(`Function.surjInv`)。原典は `U_K ⊆ Γ_K^ab` と
   書いて `Γ_K` そのものへは降ろさない。我々は `IsUniformizing` の型
   (`ρ : K.absGal →* Eˣ`)に合わせるために `Γ_K` へ持ち上げた。
   合成 `artinUnitChar ∘ artinToGal` は選択に依らない(`artinUnitChar_artinToGal`)。
2. **`Art_π` は素元 `π` と Lubin-Tate 級数 `f` の選択に依存する**
   (`Found/PGC/ArtinMap.lean` の射程)。ゆえに実物も 1 つには決まらず、
   `ArtinDatum K` という**構造体**にして「どの Artin データでも」と量化する形にした。
   ★`Classical.arbitrary` で 1 つ選んだ版(`realToGal`)も置くが、
   固定された形の主張(`Cor33Pinned`)は `ArtinDatum` を量化する方を採る
   ——選択に依存する対象を `Classical.arbitrary` で固定すると、
   `α` で結ばれた 2 つの体で選択が対応する保証が無いからである。
3. **Hodge-Tate の定義を Tate 捻りではなく固有空間で書いた**(上記)。
   数学的には同値だが、原典 [1] Serre III §1.2 の字面とは異なる。
4. **`∑_i d_V(i) = dim V` を `∃ S : Finset ℤ, ∑_{i ∈ S} d_V(i) = dim V` と書いた。**
   `d_V(i) ≥ 0` かつ総和 `≤ dim V` なので同値だが、後者の 2 つは証明していない。
5. `Representation ℚ_[p] Γ_K V` を使った(原典の「continuous Q_p[Γ_K]-vector space」の
   **連続性は落としている**)。`Skeleton` 側の `[SMul K.absGal V]` より強い。

## 触っていないもの

`Skeleton/` には 1 行も触っていない。配線(`Skeleton` の statement をどう直すか)は
人の判断待ちである。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC ABC3.Interface.PGC
open scoped NormedField Valued Classical TensorProduct

/-! ## §1 抽象核

★以下 2 本には分岐・付値・Galois・p 進の語彙が 1 つも出てこない。 -/

section AbstractCore

variable {G H M A : Type*}

/-- 抽象核 —— 全射 `r` の切断で持ち上げてから `r` を掛け直せば元に戻る。 -/
theorem apply_restrict_surjInv {r : G → H} (hr : Function.Surjective r) (a : M → H)
    (ψ : H → A) (m : M) : ψ (r (Function.surjInv hr (a m))) = ψ (a m) :=
  congrArg ψ (Function.surjInv_eq hr (a m))

/-- 抽象核 —— `ψ ∘ t` が単射写像 `c` に一致すれば `t` は単射。 -/
theorem injective_of_comp_eq {N : Type*} {t : M → G} {ψ : G → N} {c : M → N}
    (hc : Function.Injective c) (h : ∀ m, ψ (t m) = c m) : Function.Injective t :=
  fun a b hab => hc (by rw [← h a, ← h b, hab])

/-- 抽象核 —— 単射写像は、始域に相異なる 2 点があれば定数写像でない。 -/
theorem ne_const_of_injective {t : M → G} (ht : Function.Injective t) {a b : M} (hab : a ≠ b)
    (g : G) : t ≠ fun _ => g :=
  fun h => hab (ht (by rw [h]))

end AbstractCore

variable {p : ℕ} [Fact p.Prime]

/-! ## §2 `U_K = {x : ‖x‖ = 1}` -/

/-- `‖x‖ = 1` なら `x` は `𝒪_K` の単数。 -/
noncomputable def normOneToUnits (K : PAdicLocalField p)
    (x : {y : K.carrier // ‖y‖ = (1 : ℝ)}) : (𝒪[K.carrier])ˣ := by
  have hmem : (x : K.carrier) ∈ 𝒪[K.carrier] := by
    rw [Valuation.mem_integer_iff]
    have hv : Valued.v (x : K.carrier) = (‖(x : K.carrier)‖₊ : NNReal) := NNReal.eq rfl
    rw [hv]
    have h1 : ‖(x : K.carrier)‖ ≤ 1 := le_of_eq x.2
    exact_mod_cast h1
  refine (?_ : IsUnit (⟨(x : K.carrier), hmem⟩ : 𝒪[K.carrier])).unit
  rw [Valued.integer.isUnit_iff_norm_eq_one]
  exact x.2

@[simp] theorem coe_normOneToUnits (K : PAdicLocalField p)
    (x : {y : K.carrier // ‖y‖ = (1 : ℝ)}) :
    ((normOneToUnits K x : 𝒪[K.carrier]) : K.carrier) = (x : K.carrier) := rfl

/-- `U_K = {‖x‖ = 1}` は開(超距離)。 -/
theorem isOpen_normEqOne (K : PAdicLocalField p) : IsOpen {x : K.carrier | ‖x‖ = 1} := by
  rw [Metric.isOpen_iff]
  intro x hx
  refine ⟨1, one_pos, ?_⟩
  intro y hy
  simp only [Metric.mem_ball, dist_eq_norm] at hy
  simp only [Set.mem_setOf_eq] at hx ⊢
  have hne : ‖y - x‖ ≠ ‖x‖ := by rw [hx]; exact ne_of_lt hy
  have hyx : y = (y - x) + x := by ring
  rw [hyx, IsUltrametricDist.norm_add_eq_max_of_norm_ne_norm (S := K.carrier) hne, hx]
  exact max_eq_right (le_of_lt (by rw [← hx] at hy ⊢; exact hy))

/-- `U_K` には `1` 以外の元がある(位相が非離散だから)。 -/
theorem exists_normEqOne_ne_one (K : PAdicLocalField p) :
    ∃ x : K.carrier, ‖x‖ = 1 ∧ x ≠ 1 := by
  obtain ⟨δ, hδpos, hδlt⟩ := NormedField.exists_norm_lt K.carrier (one_pos (α := ℝ))
  have hne : ‖(1 : K.carrier)‖ ≠ ‖δ‖ := by rw [norm_one]; exact ne_of_gt hδlt
  refine ⟨1 + δ, ?_, ?_⟩
  · rw [IsUltrametricDist.norm_add_eq_max_of_norm_ne_norm (S := K.carrier) hne, norm_one]
    exact max_eq_left (le_of_lt hδlt)
  · intro h
    have : δ = 0 := by linear_combination h
    rw [this, norm_zero] at hδpos
    exact lt_irrefl 0 hδpos

/-! ## §3 `Γ_K ↠ Gal(K^ab/K)` -/

/-- ★`Γ_K → Gal(K^ab/K)` は全射(`K^ab/K` も `K̄/K` も正規だから)。 -/
theorem restrictNormalHom_abelianClosure_surjective (K : PAdicLocalField p) :
    Function.Surjective (AlgEquiv.restrictNormalHom
      ((abelianClosure K : IntermediateField K.carrier K.closure) : Type _) :
        K.absGal →* (abelianClosure K ≃ₐ[K.carrier] abelianClosure K)) :=
  AlgEquiv.restrictNormalHom_surjective (F := K.carrier) K.closure

/-! ## §4 Artin データ -/

/-- `Found/PGC/ArtinMap.lean::exists_artinMap` が与える組を構造体にしたもの。

★`Art_π` は素元 `π` と Lubin-Tate 級数 `f` の選択に依存するので、
「1 つの実物」ではなく「実物の集まり」として扱う。 -/
structure ArtinDatum (K : PAdicLocalField p) where
  /-- `Gal(K^ab/K) ≅ 𝒪_K^× × Gal(K^ur/K)` -/
  Φ : (abelianClosure K ≃ₐ[K.carrier] abelianClosure K) ≃* ((𝒪[K.carrier])ˣ × unramGal K)
  /-- `Art_π : K^× →* Gal(K^ab/K)` -/
  Art : (K.carrier)ˣ →* (abelianClosure K ≃ₐ[K.carrier] abelianClosure K)
  /-- `𝒪_K^×` の上での値(原典 Yoshida08 Proposition 4.7(ii) の `j = 0` の場合) -/
  art_unit : ∀ u : (𝒪[K.carrier])ˣ, Art (unitsToCarrier K u) = Φ.symm (u, 1)

def ArtinDatum.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 9, item := "Definition 4.10", sectionId := "def-4-10" }

instance nonempty_artinDatum (K : PAdicLocalField p) : Nonempty (ArtinDatum K) := by
  obtain ⟨Φ, Art, ϖ, -, hunit, -⟩ := exists_artinMap K
  exact ⟨⟨Φ, Art, hunit⟩⟩

/-- ★★実物の `toGal` —— `U_K → K^× --Art--> Gal(K^ab/K)` を `Γ_K` へ持ち上げたもの。 -/
noncomputable def artinToGal (K : PAdicLocalField p) (D : ArtinDatum K) :
    {y : K.carrier // ‖y‖ = (1 : ℝ)} → K.absGal :=
  fun x => Function.surjInv (restrictNormalHom_abelianClosure_surjective K)
    (D.Art (unitsToCarrier K (normOneToUnits K x)))

def artinToGal.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Definition 3.2", sectionId := "def-3-2" }

/-- ★実物の指標 `Γ_K ↠ Gal(K^ab/K) ≅ 𝒪_K^× × Gal(K^ur/K) → 𝒪_K^× ↪ K^×`。 -/
noncomputable def artinUnitChar (K : PAdicLocalField p) (D : ArtinDatum K) :
    K.absGal →* (K.carrier)ˣ :=
  (unitsToCarrier K).comp ((MonoidHom.fst _ _).comp
    (D.Φ.toMonoidHom.comp (AlgEquiv.restrictNormalHom
      ((abelianClosure K : IntermediateField K.carrier K.closure) : Type _))))

/-- ★★持ち上げの選択に依らないこと —— `artinUnitChar ∘ artinToGal` は `U_K ↪ K^×`。 -/
theorem artinUnitChar_artinToGal (K : PAdicLocalField p) (D : ArtinDatum K)
    (x : {y : K.carrier // ‖y‖ = (1 : ℝ)}) :
    artinUnitChar K D (artinToGal K D x) = unitsToCarrier K (normOneToUnits K x) := by
  have h : (AlgEquiv.restrictNormalHom
      ((abelianClosure K : IntermediateField K.carrier K.closure) : Type _) :
        K.absGal →* (abelianClosure K ≃ₐ[K.carrier] abelianClosure K)) (artinToGal K D x)
      = D.Art (unitsToCarrier K (normOneToUnits K x)) :=
    Function.surjInv_eq (restrictNormalHom_abelianClosure_surjective K) _
  simp only [artinUnitChar, MonoidHom.comp_apply, h, D.art_unit, MulEquiv.coe_toMonoidHom,
    MulEquiv.apply_symm_apply]
  rfl

@[simp] theorem coe_artinUnitChar_artinToGal (K : PAdicLocalField p) (D : ArtinDatum K)
    (x : {y : K.carrier // ‖y‖ = (1 : ℝ)}) :
    ((artinUnitChar K D (artinToGal K D x) : (K.carrier)ˣ) : K.carrier) = (x : K.carrier) := by
  rw [artinUnitChar_artinToGal, unitsToCarrier_val, coe_normOneToUnits]

theorem artinToGal_injective (K : PAdicLocalField p) (D : ArtinDatum K) :
    Function.Injective (artinToGal K D) :=
  injective_of_comp_eq (c := fun x : {y : K.carrier // ‖y‖ = (1 : ℝ)} => (x : K.carrier))
    (ψ := fun g => ((artinUnitChar K D g : (K.carrier)ˣ) : K.carrier))
    Subtype.val_injective (coe_artinUnitChar_artinToGal K D)

/-- ★★実物の `toGal` は定数写像ではない。

★`Check/PGC/FreeTermFunctionRefutation.lean::badToGal` は捻り側で定数 `1` を返すもの
だったから、`toGal` を実物に固定すればあの反例は作れない。 -/
theorem artinToGal_ne_const (K : PAdicLocalField p) (D : ArtinDatum K) (g : K.absGal) :
    artinToGal K D ≠ fun _ => g := by
  obtain ⟨x, hx, hxne⟩ := exists_normEqOne_ne_one K
  refine ne_const_of_injective (artinToGal_injective K D)
    (a := (⟨x, hx⟩ : {y : K.carrier // ‖y‖ = (1 : ℝ)}))
    (b := (⟨1, norm_one⟩ : {y : K.carrier // ‖y‖ = (1 : ℝ)})) ?_ g
  exact fun h => hxne (congrArg Subtype.val h)

/-! ## §5 実物の `toGal` で Definition 3.2 が真になる -/

/-- ★★★★★実物の Artin 写像と実物の指標について、Definition 3.2 が**真になる**。

開部分群は `I = U_K` そのもの、体準同型は `ι = id : K → K`、`E = K`。
★仮定ゼロ(`ArtinDatum K` は `nonempty_artinDatum` で常に非空)。 -/
theorem isUniformizing_artin (K : PAdicLocalField p) (D : ArtinDatum K) :
    IsUniformizing K K.carrier (artinToGal K D) (artinUnitChar K D) := by
  refine ⟨{x : K.carrier | ‖x‖ = 1}, subset_rfl, isOpen_normEqOne K, by simp, ?_, ?_,
    RingHom.id K.carrier, ?_⟩
  · intro a ha b hb
    simp only [Set.mem_setOf_eq] at *
    rw [norm_mul, ha, hb, one_mul]
  · intro a ha
    simp only [Set.mem_setOf_eq] at *
    rw [norm_inv, ha, inv_one]
  · intro x hx
    exact coe_artinUnitChar_artinToGal K D ⟨x, hx⟩

def isUniformizing_artin.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Definition 3.2", sectionId := "def-3-2" }

#print axioms isUniformizing_artin
#print axioms artinToGal_ne_const

/-! ### 選択で 1 つ固定した版

★`Cor33Pinned` はこちらではなく `ArtinDatum` を量化する方を使う(逸脱の記録 2)。 -/

/-- 選択公理で Artin データを 1 つ選ぶ。 -/
noncomputable def realArtinDatum (K : PAdicLocalField p) : ArtinDatum K := Classical.arbitrary _

/-- `toGal` を `K` だけから決まるものに固定した版。 -/
noncomputable def realToGal (K : PAdicLocalField p) :
    {y : K.carrier // ‖y‖ = (1 : ℝ)} → K.absGal :=
  artinToGal K (realArtinDatum K)

/-- 対応する指標。 -/
noncomputable def realUnitChar (K : PAdicLocalField p) : K.absGal →* (K.carrier)ˣ :=
  artinUnitChar K (realArtinDatum K)

theorem isUniformizing_real (K : PAdicLocalField p) :
    IsUniformizing K K.carrier (realToGal K) (realUnitChar K) :=
  isUniformizing_artin K (realArtinDatum K)

/-! ## §6 Corollary 3.3 を実物に固定した形 -/

/-- ★★[pGC] Corollary 3.3 の、`toGal` を実物に固定した形。

原文 (pGC p.6):
> Given a continuous E[Γ_K]-module V of E-dimension 1, the issue of whether or not V is
> uniformizing can be determined entirely group-theoretically from the filtered group Γ_K.

`Skeleton/PGC/Section3.lean::cor_3_3` との違いは 1 点だけ ——
`toGal` が自由なパラメータではなく `artinToGal`(Artin 写像の持ち上げ)であること。
★これは `def` であって主張の証明ではない。証明は未着手。 -/
def Cor33Pinned (RF : RamificationFiltration p) (E : Type) [Field E] [Algebra ℚ_[p] E] : Prop :=
  ∀ {K K' : PAdicLocalField p}
    (α : FilteredGroup.Iso (filteredGroupOf RF K) (filteredGroupOf RF K'))
    (DK : ArtinDatum K) (DK' : ArtinDatum K')
    (ρ : K.absGal →* Eˣ) (ρ' : K'.absGal →* Eˣ)
    (_hρ : ∀ g : K.absGal, ρ' (α.equiv g) = ρ g),
    IsUniformizing K E (artinToGal K DK) ρ ↔ IsUniformizing K' E (artinToGal K' DK') ρ'

def Cor33Pinned.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.3", sectionId := "cor-3-3" }

/-! ## §7 `ℂ_K` を `ℚ_p`-代数として見る -/

/-- `K^al` は `ℚ_p` 上のノルム代数。

★既存の `Algebra ℚ_[p] K.closure` を再利用して `norm_smul_le` だけ足しているので
菱形は作っていない。これが無いと `Module ℚ_[p] (CompKbar K)` が合成できない
(モジュール docstring §2 の実測)。 -/
noncomputable scoped instance closureNormedAlgebraQp (K : PAdicLocalField p) :
    NormedAlgebra ℚ_[p] K.closure :=
  { (inferInstance : Algebra ℚ_[p] K.closure) with
    norm_smul_le := fun r x => by
      rw [Algebra.smul_def, norm_mul, IsScalarTower.algebraMap_apply ℚ_[p] K.carrier K.closure,
        norm_algebraMap_closure, ABC3.Found.PGC.norm_algebraMap] }

/-- `Γ_K` は `ℂ_K` の `ℚ_p`-元を固定する。 -/
theorem smul_algebraMap_qp (K : PAdicLocalField p) (σ : K.absGal) (c : ℚ_[p]) :
    σ • (algebraMap ℚ_[p] (CompKbar K) c) = algebraMap ℚ_[p] (CompKbar K) c := by
  rw [IsScalarTower.algebraMap_apply ℚ_[p] K.carrier (CompKbar K), smul_algebraMap]

/-- `Γ_K` の作用は `ℚ_p`-スカラーと可換。 -/
theorem smul_qp_smul (K : PAdicLocalField p) (σ : K.absGal) (c : ℚ_[p]) (z : CompKbar K) :
    σ • (c • z) = c • (σ • z) := by
  rw [Algebra.smul_def, Algebra.smul_def, smul_mul', smul_algebraMap_qp]

/-! ## §8 `Γ_K` の `ℂ_K` への `ℚ_p`-線型表現 -/

/-- `Γ_K → GL_{ℚ_p}(ℂ_K)`。 -/
noncomputable def compRep (K : PAdicLocalField p) :
    Representation ℚ_[p] K.absGal (CompKbar K) where
  toFun σ :=
    { toFun := fun z => σ • z
      map_add' := fun z w => smul_add σ z w
      map_smul' := fun c z => smul_qp_smul K σ c z }
  map_one' := by ext z; exact one_smul _ z
  map_mul' σ τ := by ext z; exact mul_smul σ τ z

@[simp] theorem compRep_apply (K : PAdicLocalField p) (σ : K.absGal) (z : CompKbar K) :
    compRep K σ z = σ • z := rfl

/-! ## §9 円分指標 -/

/-- `χ_K : Γ_K →* ℤ_[p]ˣ`(mathlib の `cyclotomicCharacter` を `Γ_K` の上に載せたもの)。 -/
noncomputable def cycloChar (K : PAdicLocalField p) : K.absGal →* ℤ_[p]ˣ :=
  (cyclotomicCharacter K.closure p).comp
    { toFun := AlgEquiv.toRingEquiv, map_one' := rfl, map_mul' := fun _ _ => rfl }

/-- `χ_K(σ)^i` を `ℚ_p` のスカラーとして見たもの。 -/
noncomputable def cycloScalar (K : PAdicLocalField p) (σ : K.absGal) (i : ℤ) : ℚ_[p] :=
  ((cycloChar K σ ^ i : ℤ_[p]ˣ) : ℤ_[p])

@[simp] theorem cycloScalar_zero (K : PAdicLocalField p) (σ : K.absGal) :
    cycloScalar K σ 0 = 1 := by
  simp [cycloScalar]

/-! ## §10 `d_V(i)` と Hodge-Tate 性 -/

variable (K : PAdicLocalField p) (V : Type) [AddCommGroup V] [Module ℚ_[p] V]

/-- 対角作用は `K`-スカラーと可換(`Γ_K` が `K` を固定するから)。 -/
theorem tprod_smul_comm (ρ : Representation ℚ_[p] K.absGal V) (σ : K.absGal) (c : K.carrier)
    (z : CompKbar K ⊗[ℚ_[p]] V) :
    (Representation.tprod (compRep K) ρ) σ (c • z)
      = c • (Representation.tprod (compRep K) ρ) σ z := by
  induction z using TensorProduct.induction_on with
  | zero => simp
  | tmul x v =>
      rw [TensorProduct.smul_tmul']
      simp only [Representation.tprod_apply, TensorProduct.map_tmul, compRep_apply]
      rw [TensorProduct.smul_tmul', smul_comm]
  | add a b ha hb => rw [smul_add, map_add, map_add, ha, hb, smul_add]

/-- ★★`ℂ_K ⊗_{ℚ_p} V` の `χ^i`-固有空間 —— 原典の `d_V(i)` が数える `K`-ベクトル空間。

`V(−i) ⊗ ℂ_K` の `Γ_K`-不変元と同じもの(逸脱の記録 3)。 -/
def hodgeTateWeightSpace (ρ : Representation ℚ_[p] K.absGal V) (i : ℤ) :
    Submodule K.carrier (CompKbar K ⊗[ℚ_[p]] V) where
  carrier := {z | ∀ σ : K.absGal,
    (Representation.tprod (compRep K) ρ) σ z = cycloScalar K σ i • z}
  add_mem' {a b} ha hb σ := by rw [map_add, ha σ, hb σ, smul_add]
  zero_mem' σ := by rw [map_zero, smul_zero]
  smul_mem' c z hz σ := by rw [tprod_smul_comm K V ρ σ c z, hz σ, smul_comm]

def hodgeTateWeightSpace.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

theorem mem_hodgeTateWeightSpace (ρ : Representation ℚ_[p] K.absGal V) (i : ℤ)
    (z : CompKbar K ⊗[ℚ_[p]] V) :
    z ∈ hodgeTateWeightSpace K V ρ i ↔
      ∀ σ : K.absGal, (Representation.tprod (compRep K) ρ) σ z = cycloScalar K σ i • z :=
  Iff.rfl

/-- ★★原典の不変量 `d_V(i)` —— `K` 上の次元。 -/
noncomputable def hodgeTateDim (ρ : Representation ℚ_[p] K.absGal V) (i : ℤ) : ℕ :=
  Module.finrank K.carrier (hodgeTateWeightSpace K V ρ i)

def hodgeTateDim.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★`V` が Hodge-Tate であること —— `∑_i d_V(i) = dim_{ℚ_p} V`(逸脱の記録 4)。 -/
def IsHodgeTate (ρ : Representation ℚ_[p] K.absGal V) : Prop :=
  ∃ S : Finset ℤ, (∑ i ∈ S, hodgeTateDim K V ρ i) = Module.finrank ℚ_[p] V

def IsHodgeTate.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ### 退化していないことの検査 -/

/-- 自明表現の重み `0` の空間には `1 ⊗ 1` が入っている。 -/
theorem one_tmul_one_mem_weightSpace_zero :
    (1 : CompKbar K) ⊗ₜ[ℚ_[p]] (1 : ℚ_[p])
      ∈ hodgeTateWeightSpace K ℚ_[p] (Representation.trivial ℚ_[p] K.absGal ℚ_[p]) 0 := by
  intro σ
  simp only [Representation.tprod_apply, TensorProduct.map_tmul, compRep_apply,
    Representation.trivial_apply, cycloScalar_zero, one_smul]
  rw [smul_one]

/-- ★重み空間は `0` でない(定義が空虚でないことの確認)。

★★これ以上は言えない —— `d_V(0) = 1` には `ℂ_K^{Γ_K} = K`(Ax–Sen–Tate)が要る。 -/
theorem weightSpace_zero_ne_bot :
    hodgeTateWeightSpace K ℚ_[p] (Representation.trivial ℚ_[p] K.absGal ℚ_[p]) 0 ≠ ⊥ := by
  intro h
  have hz : (1 : CompKbar K) ⊗ₜ[ℚ_[p]] (1 : ℚ_[p]) = 0 :=
    (Submodule.mem_bot (R := K.carrier)).mp (h ▸ one_tmul_one_mem_weightSpace_zero K)
  have h2 := congrArg (TensorProduct.rid ℚ_[p] (CompKbar K)) hz
  rw [TensorProduct.rid_tmul, map_zero, one_smul] at h2
  exact one_ne_zero h2

#print axioms weightSpace_zero_ne_bot

/-! ## §11 Corollary 3.1 を実物に固定した形 -/

/-- ★★[pGC] Corollary 3.1 の、`isHodgeTate` を実物に固定した形。

原文 (pGC p.6):
> Given a continuous Q_p[Γ_K]-vector space of finite Q_p-dimension, the issue of whether or
> not V is Hodge-Tate (as well as the invariants d_V(i)) can be determined entirely
> group-theoretically from the filtered group Γ_K.

`Skeleton/PGC/Section3.lean::cor_3_1` との違いは 1 点だけ ——
`isHodgeTate` が自由な述語ではなく `IsHodgeTate`(`CompKbar` から作った実物)であること。
★これは `def` であって主張の証明ではない。証明は未着手。 -/
def Cor31Pinned (RF : RamificationFiltration p) : Prop :=
  ∀ {K K' : PAdicLocalField p}
    (α : FilteredGroup.Iso (filteredGroupOf RF K) (filteredGroupOf RF K'))
    (W : Type) [AddCommGroup W] [Module ℚ_[p] W] [FiniteDimensional ℚ_[p] W]
    (ρ : Representation ℚ_[p] K.absGal W) (ρ' : Representation ℚ_[p] K'.absGal W)
    (_hcompat : ∀ g : K.absGal, ρ' (α.equiv g) = ρ g),
    IsHodgeTate K W ρ ↔ IsHodgeTate K' W ρ'

def Cor31Pinned.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end ABC3.Found.PGC
