import ABC3.Found.PGC.UnramifiedCompletion

/-!
# `𝓀_{K^ur} = 𝔽̄_q` と算術 Frobenius(`sorry` 無し)

経路 Λ の節点 Λ5b。Λ5(`UnramifiedZhat.lean` / `UnramifiedCompletion.lean`)は
**意図的に剰余体へ降りずに** `Gal(K^ur/K) ≃ₜ* Ẑ` と `K̂^{ur}` を作った。
Λ6(Dwork の補題 = `Art_π` の π 非依存性)はそこを通れない——
逐次近似の各段で

* `𝓀_{K^ur}` が**代数閉体**であること(`x^{q-1} = u` と `c^q - c = a` の可解性)、
* Frobenius が剰余体上で **`x ↦ x^q`** であること

を使うからである。本ファイルはその 3 つを埋める。

## 到達点

| | 主張 |
|---|---|
| 目標 1 | `card_residueField_unramLevel`:`K^ur` の次数 `N` の段の剰余体は `q^N` 元 |
| 目標 2 | `isAlgClosed_unramResidueField`:`𝓀_{K^ur}` は代数閉体(= `𝔽̄_q`) |
| 目標 3 | `exists_arithFrobenius`:`Gal(K^ur/K)` に剰余体上 `x ↦ x^q` の元がある |
| Λ6 の入口 | `isAlgClosed_residueField_unramifiedCompletionInt` / `exists_arithFrobenius_completion` / `exists_pow_eq_completion` / `exists_pow_card_sub_self_eq_completion` |

## ★★★訂正:`coherentFrobenius` は算術 Frobenius **ではない**

`UnramifiedZhat.lean` の `exists_coherentFrobenius` は
`unramLevelGeneratorSet K N`(= `G/P_N ≅ ℤ/N` を生成する元)の共通部分から
コンパクト性で 1 点を取るので、得られるのは `Ẑ ≅ Gal(K^ur/K)` の
**位相的生成元一般**である。位相的生成元の全体は Frobenius の `Ẑ^×`-軌道
なので、`coherentFrobenius K` の剰余体上の作用は一般に `x ↦ x^{q^k}`
(`k ∈ Ẑ^×`)であって `x ↦ x^q` とは限らない。
したがって「`coherentFrobenius` が剰余体上 `x ↦ x^q`」という形の主張は
**偽になりうる**。本ファイルは代わりに**算術 Frobenius の存在**
(`exists_arithFrobenius`)を、剰余体の条件を直接課したコンパクト性で示す。
★Dwork の補題は算術 Frobenius でなければ回らない(`σ(ξ)/ξ = u` の剰余方程式
`ξ̄^{q^k}/ξ̄ = ū` は `k ∈ Ẑ` では多項式方程式にならない)。

## 筋(剰余体を `K̄` 側で捕まえる)

古典的には「`𝓀_{K^ur}` は `𝔽_{q^n}` の合併だから代数閉」と言うが、
それは有限部分体の有向系を組み、多項式の係数を 1 つの段へ持ち上げる
配管を要する。本ファイルは**代数閉包の剰余体を経由する**近道を取る。

1. **`𝓀_{K̄}` は代数閉**(`isAlgClosed_residueField_closureInt`)。
   モニック多項式を `𝒪_{K̄}` へ持ち上げ、`K̄` の中で根を取り、
   「モニックで係数が整な多項式の根は整」(`norm_le_one_of_monic_eval_eq_zero`、
   超距離不等式だけの初等的な補題)から根が `𝒪_{K̄}` に入ることを見て還元する。
2. **`𝓀_{K^ur} ≅ 𝓀_{K̄}`**(`unramResidueFieldEquivClosure`)。単射は体からの
   環準同型だから。全射は**根の数え上げ**:`t ∈ 𝓀_{K̄}` は `t^{q^n} = t` を満たし
   (`K(t)` の剰余体が有限だから)、`𝓀_{K^ur}` の中の `q^n` 元の部分体
   (目標 1)は `X^{q^n} - X` の根を**すべて**尽くす(根はたかだか `q^n` 個)。
3. **算術 Frobenius**。`Γ_K` の中で
   `B_x := {g | g は 𝒪_{K̄} ∩ K(x) の剰余に q 乗として作用する}` を考える。
   `B_x` は開部分群 `K(x).fixingSubgroup` による右移動で閉じるので **clopen**、
   `exists_frobenius`(有限次不分岐拡大の Frobenius)と
   `AlgEquiv.restrictNormalHom_surjective` で**空でない**、次数の整除で**有向**。
   `Γ_K` はコンパクトなので共通部分に元 `g` があり、それを `K^ur` へ制限する。

## mathlib の在庫でどこまで済んだか(2026-09-06 の実測)

| 要るもの | 在庫 |
|---|---|
| モニック既約多項式の持ち上げ | `Polynomial.lifts_and_natDegree_eq_and_monic` |
| 代数閉性の判定 | `IsAlgClosed.of_exists_root` / `of_ringEquiv` |
| 有限体の `x^{|F|} = x` | `FiniteField.pow_card` |
| 根の個数 | `Polynomial.card_roots'` + `Finset.eq_of_subset_of_card_le` |
| 超距離の有限和 | `IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg` |
| Galois の制限が全射 | `AlgEquiv.restrictNormalHom_surjective` |
| `Γ_K` のコンパクト性 | `InfiniteGalois`(`IsGalois K.carrier K.closure` から `infer_instance`) |
| 有限次の Frobenius | `UnramifiedExtension.exists_frobenius`(本プロジェクト) |
| 完備化は剰余体を変えない | `UnramifiedCompletion.residueFieldEquiv`(本プロジェクト) |

★本ファイルは docstring 込みで約 1,015 行(証明本体は約 750 行)。
目標 1・2・3 の 3 つを 1 ファイルで埋めた分、行数そのものは前回までの
節点より大きい。ただし**当初の見立て(有限部分体の有向系を組み、多項式の
係数を 1 つの段へ持ち上げる)は取らずに済んだ**——「代数閉包の剰余体を
経由する」という筋の取り替えが効いており、削れたのは mathlib の在庫よりも
**設計の選択**による。

## ★設計上の注意(守ったこと)

* **`inertia` を経由していない**。不分岐側は `unramifiedClosure K` という体として
  直接扱っており、`Interface` の `SubgroupCorrespondence` / `ResidueCardinality` は
  本ファイルの主張にも証明にも現れない(Corollary 1.3 との循環を避けるため)。
* **`Abelianization` を使っていない**。
* **結論に自由なパラメータを出していない**——`exists_arithFrobenius` の型は
  `K` にしか依存せず、`σ` は `∃` の内側にある。
* **`sorry` 本体の `def` を作っていない**(D21)。
* 今日入った 9 ファイルは**変更していない**(import のみ)。

## 逸脱(記録)

* 「`𝓀_{K^ur} = 𝔽̄_q`」を、有限部分体の合併としてではなく
  **`𝓀_{K̄}` との同型**として与えた(上の筋 2)。数学的には同値だが、
  `𝔽_{q^n}` の有向系や `GaloisField p n` との同型は**作っていない**。
  下流が「`𝔽_{q^n}` そのもの」を要求したら節点を立てること。
* `𝒪_{K̄}`(`closureInt`)は `Valued` を経由せずノルム `≤ 1` の部分環として
  作った。`K.closure` に `Valued` インスタンスを足すと、他所の
  `NormedField.toValued` と菱形になる恐れがあるため
  (`UnramifiedCompletion.lean` と同じ判断)。
* 目標 3 の `σ` が**位相的生成元でもある**ことは主張していない。
  実際そうであるが(剰余体の Frobenius は各段で位数 `N`)、Dwork には
  不要なので示していない。必要になったら節点を立てること。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped Valued

variable {p : ℕ} [Fact p.Prime]

/-! ## 1. 超距離ノルム体の単位球

`adjoinIntegers` も `closureInt` も「ノルム `≤ 1` の部分環」なので、
局所環であることと極大イデアルの記述は一度だけ書いておけばよい。 -/

section UnitBall

variable {F : Type*} [NormedField F] {S : Subring F}

/-- ノルムが `1` の整数は単元(逆元のノルムも `1`)。 -/
theorem isUnit_of_norm_eq_one' (hS : ∀ y : F, y ∈ S ↔ ‖y‖ ≤ 1) {w : S}
    (h : ‖(w : F)‖ = 1) : IsUnit w := by
  have hne : (w : F) ≠ 0 := by intro h0; rw [h0, norm_zero] at h; exact zero_ne_one h
  have hmem : (w : F)⁻¹ ∈ S := (hS _).mpr (by rw [norm_inv, h]; simp)
  exact ⟨⟨w, ⟨(w : F)⁻¹, hmem⟩, Subtype.ext (mul_inv_cancel₀ hne),
    Subtype.ext (inv_mul_cancel₀ hne)⟩, rfl⟩

/-- 単元のノルムは `1`(自分と逆元のノルムが両方 `≤ 1` で積が `1`)。 -/
theorem norm_eq_one_of_isUnit' (hS : ∀ y : F, y ∈ S ↔ ‖y‖ ≤ 1) {w : S}
    (h : IsUnit w) : ‖(w : F)‖ = 1 := by
  obtain ⟨u, rfl⟩ := h
  have h1 : ((u : S) : F) * ((↑u⁻¹ : S) : F) = 1 := by
    rw [← Subring.coe_mul, ← Units.val_mul, mul_inv_cancel, Units.val_one, Subring.coe_one]
  have h2 : ‖((u : S) : F)‖ * ‖((↑u⁻¹ : S) : F)‖ = 1 := by rw [← norm_mul, h1, norm_one]
  have ha : ‖((u : S) : F)‖ ≤ 1 := (hS _).mp (u : S).2
  have hb : ‖((↑u⁻¹ : S) : F)‖ ≤ 1 := (hS _).mp (↑u⁻¹ : S).2
  nlinarith [norm_nonneg (((u : S) : F)), norm_nonneg ((↑u⁻¹ : S) : F)]

theorem isUnit_iff_norm_eq_one' (hS : ∀ y : F, y ∈ S ↔ ‖y‖ ≤ 1) (w : S) :
    IsUnit w ↔ ‖(w : F)‖ = 1 :=
  ⟨norm_eq_one_of_isUnit' hS, isUnit_of_norm_eq_one' hS⟩

/-- **超距離ノルム体の単位球は局所環**——非単元(ノルム `< 1`)は
非アルキメデス三角不等式で加法について閉じる。

退化の自己検査:`IsUltrametricDist` を落とすと**偽**(ℝ の単位球は局所環でない)。 -/
theorem isLocalRing_of_norm_le_one [IsUltrametricDist F] (hS : ∀ y : F, y ∈ S ↔ ‖y‖ ≤ 1) :
    IsLocalRing S := by
  haveI : Nontrivial S := ⟨⟨0, 1, by simp⟩⟩
  refine IsLocalRing.of_nonunits_add ?_
  intro a b ha hb
  rw [mem_nonunits_iff, isUnit_iff_norm_eq_one' hS] at *
  have hab : ‖((a + b : S) : F)‖ ≤ max ‖(a : F)‖ ‖(b : F)‖ := by
    rw [Subring.coe_add]; exact IsUltrametricDist.norm_add_le_max _ _
  have ha1 : ‖(a : F)‖ < 1 := lt_of_le_of_ne ((hS _).mp a.2) ha
  have hb1 : ‖(b : F)‖ < 1 := lt_of_le_of_ne ((hS _).mp b.2) hb
  exact ne_of_lt (lt_of_le_of_lt hab (max_lt ha1 hb1))

theorem mem_maximalIdeal_iff_norm_lt_one [IsUltrametricDist F] [IsLocalRing S]
    (hS : ∀ y : F, y ∈ S ↔ ‖y‖ ≤ 1) (w : S) :
    w ∈ IsLocalRing.maximalIdeal S ↔ ‖(w : F)‖ < 1 := by
  rw [IsLocalRing.mem_maximalIdeal, mem_nonunits_iff, isUnit_iff_norm_eq_one' hS]
  exact ⟨fun h => lt_of_le_of_ne ((hS _).mp w.2) h, fun h => ne_of_lt h⟩

end UnitBall

/-! ## 2. `𝒪_{K̄}` —— 代数閉包の整数環 -/

/-- `adjoinIntegers K x` は「ノルム `≤ 1`」で決まる部分環(定義そのもの)。 -/
theorem mem_adjoinIntegers_iff (K : PAdicLocalField p) (x : K.closure)
    (y : ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure))) :
    y ∈ adjoinIntegers K x ↔ ‖y‖ ≤ 1 := Iff.rfl

theorem mem_maximalIdeal_adjoinIntegers (K : PAdicLocalField p) (x : K.closure)
    (w : adjoinIntegers K x) :
    w ∈ IsLocalRing.maximalIdeal (adjoinIntegers K x)
      ↔ ‖(w : ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure)))‖ < 1 :=
  mem_maximalIdeal_iff_norm_lt_one (mem_adjoinIntegers_iff K x) w

/-- **`𝒪_{K̄}`** —— 代数閉包 `K.closure` の整数環(ノルム `≤ 1` の部分環)。

★`Valued` を経由しない。`K.closure` に `Valued` を足すと他所の
`NormedField.toValued` と菱形になるため(モジュール docstring の逸脱記録)。 -/
noncomputable def closureInt (K : PAdicLocalField p) : Subring K.closure where
  carrier := {y | ‖y‖ ≤ 1}
  mul_mem' := by
    intro a b ha hb
    simp only [Set.mem_setOf_eq] at *
    calc ‖a * b‖ = ‖a‖ * ‖b‖ := norm_mul a b
      _ ≤ 1 * 1 := mul_le_mul ha hb (norm_nonneg b) zero_le_one
      _ = 1 := by ring
  one_mem' := by simp
  add_mem' := by
    intro a b ha hb
    simp only [Set.mem_setOf_eq] at *
    calc ‖a + b‖ ≤ max ‖a‖ ‖b‖ := IsUltrametricDist.norm_add_le_max a b
      _ ≤ 1 := max_le ha hb
  zero_mem' := by simp
  neg_mem' := by
    intro a ha
    simp only [Set.mem_setOf_eq, norm_neg] at *
    exact ha

@[simp] theorem mem_closureInt (K : PAdicLocalField p) (y : K.closure) :
    y ∈ closureInt K ↔ ‖y‖ ≤ 1 := Iff.rfl

noncomputable instance isLocalRing_closureInt (K : PAdicLocalField p) :
    IsLocalRing ↥(closureInt K) :=
  isLocalRing_of_norm_le_one (mem_closureInt K)

theorem mem_maximalIdeal_closureInt (K : PAdicLocalField p) (w : ↥(closureInt K)) :
    w ∈ IsLocalRing.maximalIdeal ↥(closureInt K) ↔ ‖(w : K.closure)‖ < 1 :=
  mem_maximalIdeal_iff_norm_lt_one (mem_closureInt K) w

/-! ## 3. 整数環のあいだの環準同型

`𝒪_{K(x)} → 𝒪_{K^ur} → 𝒪_{K̄}`。★`unramifiedClosureInt` は `ValuationSubring`
なので、`(f z).val.val = z.val.val` を 1 本の `rfl` で書くと kernel が止まる
(`tools/lean-idioms.md` #59)。**`_apply` 補題を先に `rfl` で置いてから
`rw` で 1 層ずつ剥がす**。 -/

/-- `𝒪_{K(x)} → 𝒪_{K^ur}`(`K(x) ⊆ K^ur` のとき)。 -/
noncomputable def unramIntHomOfLe (K : PAdicLocalField p) {x : K.closure}
    (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure) ≤ unramifiedClosure K) :
    adjoinIntegers K x →+* ↥(unramifiedClosureInt K) where
  toFun z := ⟨⟨((z : ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure))) : K.closure),
      hle (z : ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure))).2⟩,
    (mem_unramifiedClosureInt K _).mpr z.2⟩
  map_one' := by apply Subtype.ext; apply Subtype.ext; rfl
  map_mul' _ _ := by apply Subtype.ext; apply Subtype.ext; rfl
  map_zero' := by apply Subtype.ext; apply Subtype.ext; rfl
  map_add' _ _ := by apply Subtype.ext; apply Subtype.ext; rfl

theorem unramIntHomOfLe_apply (K : PAdicLocalField p) {x : K.closure}
    (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure) ≤ unramifiedClosure K)
    (z : adjoinIntegers K x) :
    unramIntHomOfLe K hle z
      = ⟨⟨((z : ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure))) : K.closure),
          hle (z : ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure))).2⟩,
        (mem_unramifiedClosureInt K _).mpr z.2⟩ := rfl

theorem norm_unramIntHomOfLe (K : PAdicLocalField p) {x : K.closure}
    (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure) ≤ unramifiedClosure K)
    (z : adjoinIntegers K x) :
    ‖((unramIntHomOfLe K hle z : ↥(unramifiedClosureInt K)) : ↥(unramifiedClosure K))‖
      = ‖(z : ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure)))‖ := by
  rw [unramIntHomOfLe_apply]
  rfl

/-- `𝒪_{K(x)} → 𝒪_{K̄}`。 -/
noncomputable def adjToClosureInt (K : PAdicLocalField p) (x : K.closure) :
    adjoinIntegers K x →+* ↥(closureInt K) where
  toFun z := ⟨((z : ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure))) : K.closure),
    (mem_closureInt K _).mpr z.2⟩
  map_one' := by apply Subtype.ext; rfl
  map_mul' _ _ := by apply Subtype.ext; rfl
  map_zero' := by apply Subtype.ext; rfl
  map_add' _ _ := by apply Subtype.ext; rfl

theorem adjToClosureInt_apply (K : PAdicLocalField p) (x : K.closure) (z : adjoinIntegers K x) :
    adjToClosureInt K x z
      = ⟨((z : ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure))) : K.closure),
        (mem_closureInt K _).mpr z.2⟩ := rfl

@[simp] theorem adjToClosureInt_coe (K : PAdicLocalField p) (x : K.closure)
    (z : adjoinIntegers K x) :
    ((adjToClosureInt K x z : ↥(closureInt K)) : K.closure)
      = ((z : ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure))) : K.closure) := by
  rw [adjToClosureInt_apply]

theorem norm_adjToClosureInt (K : PAdicLocalField p) (x : K.closure) (z : adjoinIntegers K x) :
    ‖((adjToClosureInt K x z : ↥(closureInt K)) : K.closure)‖
      = ‖(z : ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure)))‖ := by
  rw [adjToClosureInt_apply]; rfl

/-- `𝒪_{K^ur} → 𝒪_{K̄}`。 -/
noncomputable def urToClosureInt (K : PAdicLocalField p) :
    ↥(unramifiedClosureInt K) →+* ↥(closureInt K) where
  toFun w := ⟨((w : ↥(unramifiedClosure K)) : K.closure),
    (mem_closureInt K _).mpr ((mem_unramifiedClosureInt K _).mp w.2)⟩
  map_one' := by apply Subtype.ext; rfl
  map_mul' _ _ := by apply Subtype.ext; rfl
  map_zero' := by apply Subtype.ext; rfl
  map_add' _ _ := by apply Subtype.ext; rfl

theorem urToClosureInt_apply (K : PAdicLocalField p) (w : ↥(unramifiedClosureInt K)) :
    urToClosureInt K w
      = ⟨((w : ↥(unramifiedClosure K)) : K.closure),
        (mem_closureInt K _).mpr ((mem_unramifiedClosureInt K _).mp w.2)⟩ := rfl

@[simp] theorem urToClosureInt_coe (K : PAdicLocalField p) (w : ↥(unramifiedClosureInt K)) :
    ((urToClosureInt K w : ↥(closureInt K)) : K.closure)
      = ((w : ↥(unramifiedClosure K)) : K.closure) := by
  rw [urToClosureInt_apply]

theorem norm_urToClosureInt (K : PAdicLocalField p) (w : ↥(unramifiedClosureInt K)) :
    ‖((urToClosureInt K w : ↥(closureInt K)) : K.closure)‖
      = ‖(w : ↥(unramifiedClosure K))‖ := by
  rw [urToClosureInt_apply]; rfl

instance isLocalHom_adjToClosureInt (K : PAdicLocalField p) (x : K.closure) :
    IsLocalHom (adjToClosureInt K x) := by
  constructor
  intro a ha
  rw [isUnit_iff_norm_eq_one' (mem_closureInt K), norm_adjToClosureInt] at ha
  exact isUnit_of_norm_eq_one' (mem_adjoinIntegers_iff K x) ha

instance isLocalHom_urToClosureInt (K : PAdicLocalField p) :
    IsLocalHom (urToClosureInt K) := by
  constructor
  intro a ha
  rw [isUnit_iff_norm_eq_one' (mem_closureInt K), norm_urToClosureInt] at ha
  by_contra hcon
  exact absurd ((not_isUnit_unramifiedClosureInt K a).mp hcon) (not_lt.mpr (le_of_eq ha.symm))

instance isLocalHom_unramIntHomOfLe (K : PAdicLocalField p) {x : K.closure}
    (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure) ≤ unramifiedClosure K) :
    IsLocalHom (unramIntHomOfLe K hle) := by
  constructor
  intro a ha
  by_contra hcon
  have h1 : ‖(a : ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure)))‖ < 1 :=
    lt_of_le_of_ne ((mem_adjoinIntegers_iff K x _).mp a.2)
      (fun he => hcon (isUnit_of_norm_eq_one' (mem_adjoinIntegers_iff K x) he))
  exact (not_isUnit_unramifiedClosureInt K _).mpr (by rw [norm_unramIntHomOfLe]; exact h1) ha

/-- `𝒪_{K(x)} → 𝒪_{K^ur} → 𝒪_{K̄}` は `𝒪_{K(x)} → 𝒪_{K̄}` に等しい。 -/
theorem urToClosureInt_unramIntHomOfLe (K : PAdicLocalField p) {x : K.closure}
    (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure) ≤ unramifiedClosure K)
    (z : adjoinIntegers K x) :
    urToClosureInt K (unramIntHomOfLe K hle z) = adjToClosureInt K x z := by
  rw [unramIntHomOfLe_apply, urToClosureInt_apply, adjToClosureInt_apply]

/-! ## 4. 剰余体とそのあいだの射 -/

/-- `𝓀_{K^ur}` —— 不分岐閉包の剰余体。 -/
abbrev unramResidueField (K : PAdicLocalField p) : Type :=
  IsLocalRing.ResidueField ↥(unramifiedClosureInt K)

/-- `𝓀_{K̄}` —— 代数閉包の剰余体。 -/
abbrev closureResidueField (K : PAdicLocalField p) : Type :=
  IsLocalRing.ResidueField ↥(closureInt K)

/-- `𝓀_{K^ur} → 𝓀_{K̄}`。 -/
noncomputable def unramResidueToClosure (K : PAdicLocalField p) :
    unramResidueField K →+* closureResidueField K :=
  IsLocalRing.ResidueField.map (urToClosureInt K)

/-- `𝓀_{K(x)} → 𝓀_{K̄}`。 -/
noncomputable def adjResidueToClosure (K : PAdicLocalField p) (x : K.closure) :
    IsLocalRing.ResidueField (adjoinIntegers K x) →+* closureResidueField K :=
  IsLocalRing.ResidueField.map (adjToClosureInt K x)

/-- `𝓀_{K(x)} → 𝓀_{K^ur}`(`K(x) ⊆ K^ur` のとき)。 -/
noncomputable def adjResidueToUnram (K : PAdicLocalField p) {x : K.closure}
    (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure) ≤ unramifiedClosure K) :
    IsLocalRing.ResidueField (adjoinIntegers K x) →+* unramResidueField K :=
  IsLocalRing.ResidueField.map (unramIntHomOfLe K hle)

theorem adjResidueToClosure_residue (K : PAdicLocalField p) (x : K.closure)
    (z : adjoinIntegers K x) :
    adjResidueToClosure K x (IsLocalRing.residue _ z)
      = IsLocalRing.residue _ (adjToClosureInt K x z) :=
  IsLocalRing.ResidueField.map_residue _ z

theorem adjResidueToUnram_residue (K : PAdicLocalField p) {x : K.closure}
    (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure) ≤ unramifiedClosure K)
    (z : adjoinIntegers K x) :
    adjResidueToUnram K hle (IsLocalRing.residue _ z)
      = IsLocalRing.residue _ (unramIntHomOfLe K hle z) :=
  IsLocalRing.ResidueField.map_residue _ z

theorem unramResidueToClosure_residue (K : PAdicLocalField p) (w : ↥(unramifiedClosureInt K)) :
    unramResidueToClosure K (IsLocalRing.residue _ w)
      = IsLocalRing.residue _ (urToClosureInt K w) :=
  IsLocalRing.ResidueField.map_residue _ w

/-- 三角形の可換性。 -/
theorem unramResidueToClosure_adjResidueToUnram (K : PAdicLocalField p) {x : K.closure}
    (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure) ≤ unramifiedClosure K)
    (t : IsLocalRing.ResidueField (adjoinIntegers K x)) :
    unramResidueToClosure K (adjResidueToUnram K hle t) = adjResidueToClosure K x t := by
  obtain ⟨z, rfl⟩ := IsLocalRing.residue_surjective t
  rw [adjResidueToUnram_residue, unramResidueToClosure_residue, adjResidueToClosure_residue,
    urToClosureInt_unramIntHomOfLe]

/-! ## 5. 目標 1 —— `K^ur` の次数 `N` の段の剰余体は `𝔽_{q^N}` -/

/-- `K^ur` の中の次数 `N` の不分岐中間体の生成元(`K.closure` の元として)。 -/
noncomputable def unramLevelValGen (K : PAdicLocalField p) (N : ℕ) : K.closure :=
  ((unramLevelGen K N : ↥(unramifiedClosure K)) : K.closure)

theorem isUnramified_unramLevelValGen (K : PAdicLocalField p) {N : ℕ} (hN : N ≠ 0) :
    IsUnramifiedAdjoin K (unramLevelValGen K N) := (unramLevel_spec K hN).1

theorem finrank_adjoin_unramLevelValGen (K : PAdicLocalField p) {N : ℕ} (hN : N ≠ 0) :
    Module.finrank K.carrier (IntermediateField.adjoin K.carrier
      ({unramLevelValGen K N} : Set K.closure)) = N := (unramLevel_spec K hN).2.1

/-- **★★★★★★★★★★★★★★(目標 1)`K^ur` の次数 `N` の段の剰余体は `𝔽_{q^N}`**。

`e·f = [L:K]` と不分岐性(`e = 1`)から `f = N`、有限体上の有限次元
ベクトル空間の元の個数(`residueDegree_eq_residueCard_pow`)で `q^N`。

退化の自己検査:`hN` を落とすと `unramLevel K 0` は意味を持たない
(`unramLevelGen K 0 = 0` なので `K(0) = K`、剰余体は `q` 元で `q^0 = 1` と違う)。 -/
theorem card_residueField_unramLevel (K : PAdicLocalField p) {N : ℕ} (hN : N ≠ 0) :
    Nat.card (IsLocalRing.ResidueField (adjoinIntegers K (unramLevelValGen K N)))
      = Nat.card 𝓀[K.carrier] ^ N := by
  have h1 := residueDegree_eq_residueCard_pow K (unramLevelValGen K N)
  rw [inertiaDegree_eq_finrank_of_isUnramified K _ (isUnramified_unramLevelValGen K hN),
    finrank_adjoin_unramLevelValGen K hN] at h1
  exact h1

/-! ## 6. 目標 2 の前半 —— `𝓀_{K̄}` は代数閉体 -/

open Polynomial in
/-- **モニックで係数がすべてノルム `≤ 1` の多項式の根はノルム `≤ 1`**。

`‖α‖ > 1` とすると `‖α^d‖ = ‖α‖^d` が下位項の最大値 `≤ ‖α‖^{d-1}` を超える。
超距離不等式だけを使う初等的な補題(`spectralValue` を経由しない)。 -/
theorem norm_le_one_of_monic_eval_eq_zero {F : Type*} [NormedField F] [IsUltrametricDist F]
    {P : F[X]} (hP : P.Monic) (hd : 0 < P.natDegree)
    (hc : ∀ n, ‖P.coeff n‖ ≤ 1) {α : F} (hα : P.eval α = 0) : ‖α‖ ≤ 1 := by
  by_contra hgt
  rw [not_le] at hgt
  have hpos : (0:ℝ) < ‖α‖ := lt_trans zero_lt_one hgt
  set d := P.natDegree with hdd
  have hsplit : (0:F) = α ^ d + ∑ i ∈ Finset.range d, P.coeff i * α ^ i := by
    rw [← hα, Polynomial.eval_eq_sum_range, Finset.sum_range_succ, ← hdd, hP.coeff_natDegree,
      one_mul]
    ring
  have hbound : ‖∑ i ∈ Finset.range d, P.coeff i * α ^ i‖ ≤ ‖α‖ ^ (d - 1) := by
    refine IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg (by positivity) ?_
    intro i hi
    rw [Finset.mem_range] at hi
    rw [norm_mul, norm_pow]
    calc ‖P.coeff i‖ * ‖α‖ ^ i ≤ 1 * ‖α‖ ^ i :=
          mul_le_mul_of_nonneg_right (hc i) (by positivity)
      _ = ‖α‖ ^ i := one_mul _
      _ ≤ ‖α‖ ^ (d - 1) := pow_le_pow_right₀ hgt.le (by omega)
  have heq : ‖α‖ ^ d = ‖∑ i ∈ Finset.range d, P.coeff i * α ^ i‖ := by
    rw [← norm_pow, show α ^ d = -∑ i ∈ Finset.range d, P.coeff i * α ^ i by
      rw [eq_neg_iff_add_eq_zero, ← hsplit], norm_neg]
  have hlt : ‖α‖ ^ (d - 1) < ‖α‖ ^ d := pow_lt_pow_right₀ hgt (by omega)
  linarith [heq ▸ hbound]

/-- **★★★★★★★★★★★★`𝓀_{K̄}` は代数閉体**——モニック既約多項式を `𝒪_{K̄}` へ
持ち上げ(`Polynomial.lifts_and_natDegree_eq_and_monic`)、`K̄` の中の根を取ると
(モニックで係数が整なので)その根も整で、還元するともとの多項式の根になる。

退化の自己検査:`K̄` を `K^ur` に替えるとこの証明は使えない
(`K^ur` は代数閉でないので根が取れない)——だからこそ目標 2 は
`𝓀_{K^ur} ≅ 𝓀_{K̄}` を経由する。 -/
theorem isAlgClosed_residueField_closureInt (K : PAdicLocalField p) :
    IsAlgClosed (IsLocalRing.ResidueField ↥(closureInt K)) := by
  refine IsAlgClosed.of_exists_root _ (fun f hmonic hirr => ?_)
  have hlifts : f ∈ Polynomial.lifts (IsLocalRing.residue ↥(closureInt K)) :=
    (Polynomial.lifts_iff_coeff_lifts f).mpr fun n => IsLocalRing.residue_surjective (f.coeff n)
  obtain ⟨g, hgmap, hgdeg, hgmonic⟩ := Polynomial.lifts_and_natDegree_eq_and_monic hlifts hmonic
  have hdpos : 0 < f.natDegree :=
    Polynomial.natDegree_pos_iff_degree_pos.mpr (Polynomial.degree_pos_of_irreducible hirr)
  set G : Polynomial K.closure := g.map (Subring.subtype (closureInt K)) with hG
  have hGmonic : G.Monic := hgmonic.map _
  have hGdeg : G.natDegree = f.natDegree := by rw [hG, hgmonic.natDegree_map, hgdeg]
  obtain ⟨α, hα⟩ := IsAlgClosed.exists_root G (by
    rw [Polynomial.degree_eq_natDegree hGmonic.ne_zero, hGdeg]
    exact_mod_cast hdpos.ne')
  have hαnorm : ‖α‖ ≤ 1 := by
    refine norm_le_one_of_monic_eval_eq_zero hGmonic (by rw [hGdeg]; exact hdpos) ?_ hα
    intro n
    rw [hG, Polynomial.coeff_map]
    exact (g.coeff n).2
  set a : ↥(closureInt K) := ⟨α, (mem_closureInt K α).mpr hαnorm⟩ with ha
  refine ⟨IsLocalRing.residue _ a, ?_⟩
  have hzero : g.eval a = 0 := by
    have h1 : G.eval α = ((g.eval a : ↥(closureInt K)) : K.closure) := by
      rw [hG, Polynomial.eval_map]
      exact Polynomial.eval₂_at_apply (Subring.subtype (closureInt K)) a
    rw [hα] at h1
    exact Subtype.ext h1.symm
  rw [← hgmap, Polynomial.eval_map, Polynomial.eval₂_at_apply, hzero, map_zero]

/-! ## 7. 目標 2 —— `𝓀_{K^ur} ≅ 𝓀_{K̄}` は代数閉体 -/

/-- **`𝓀_{K̄}` の各元はどれかの有限体 `𝔽_{q^n}` に入る**——`y ∈ 𝒪_{K̄}` は
有限次拡大 `K(y)` の整数環に入り、その剰余体は `q^f` 元の有限体だから。 -/
theorem exists_pow_eq_self_closureResidueField (K : PAdicLocalField p)
    (t : closureResidueField K) :
    ∃ n : ℕ, n ≠ 0 ∧ t ^ (Nat.card 𝓀[K.carrier] ^ n) = t := by
  obtain ⟨w, rfl⟩ := IsLocalRing.residue_surjective t
  set y : K.closure := (w : K.closure) with hy
  have hmem : y ∈ IntermediateField.adjoin K.carrier ({y} : Set K.closure) :=
    IntermediateField.mem_adjoin_simple_self K.carrier y
  have hnorm : (⟨y, hmem⟩ : ↥(IntermediateField.adjoin K.carrier ({y} : Set K.closure)))
      ∈ adjoinIntegers K y := (mem_closureInt K y).mp w.2
  set z : adjoinIntegers K y := ⟨⟨y, hmem⟩, hnorm⟩ with hz
  have hw : adjToClosureInt K y z = w := by
    rw [adjToClosureInt_apply]
  haveI : Fintype (IsLocalRing.ResidueField (adjoinIntegers K y)) := Fintype.ofFinite _
  have hcard : Nat.card (IsLocalRing.ResidueField (adjoinIntegers K y))
      = Nat.card 𝓀[K.carrier] ^ inertiaDegree K y := residueDegree_eq_residueCard_pow K y
  have hone : 1 < Nat.card (IsLocalRing.ResidueField (adjoinIntegers K y)) := Finite.one_lt_card
  refine ⟨inertiaDegree K y, ?_, ?_⟩
  · intro h0
    rw [h0, pow_zero] at hcard
    omega
  · have hpow : (IsLocalRing.residue (adjoinIntegers K y) z)
        ^ (Nat.card 𝓀[K.carrier] ^ inertiaDegree K y)
        = IsLocalRing.residue (adjoinIntegers K y) z := by
      rw [← hcard, Nat.card_eq_fintype_card]
      exact FiniteField.pow_card _
    calc (IsLocalRing.residue ↥(closureInt K) w) ^ (Nat.card 𝓀[K.carrier] ^ inertiaDegree K y)
        = (adjResidueToClosure K y (IsLocalRing.residue (adjoinIntegers K y) z))
            ^ (Nat.card 𝓀[K.carrier] ^ inertiaDegree K y) := by
          rw [adjResidueToClosure_residue, hw]
      _ = adjResidueToClosure K y ((IsLocalRing.residue (adjoinIntegers K y) z)
            ^ (Nat.card 𝓀[K.carrier] ^ inertiaDegree K y)) := by rw [map_pow]
      _ = adjResidueToClosure K y (IsLocalRing.residue (adjoinIntegers K y) z) := by rw [hpow]
      _ = IsLocalRing.residue ↥(closureInt K) w := by rw [adjResidueToClosure_residue, hw]

open Polynomial in
/-- **`X^m - X` の根はたかだか `m` 個**——`m` 元からなる根の集合はすべての根を尽くす。

退化の自己検査:`hTcard`(ちょうど `m` 元)を `≤ m` に弱めると**偽**。 -/
theorem mem_of_pow_eq_self_of_card_eq {F : Type*} [Field F] [DecidableEq F] {m : ℕ} (hm : 1 < m)
    {T : Finset F} (hTcard : T.card = m) (hT : ∀ s ∈ T, s ^ m = s)
    {t : F} (ht : t ^ m = t) : t ∈ T := by
  set P : F[X] := X ^ m - X with hP
  have hdeg : (X : F[X]).natDegree < (X ^ m : F[X]).natDegree := by
    rw [Polynomial.natDegree_X, Polynomial.natDegree_X_pow]; omega
  have hPmonic : P.Monic := by
    refine (Polynomial.monic_X_pow m).sub_of_left ?_
    rw [Polynomial.degree_X, Polynomial.degree_X_pow]
    exact_mod_cast hm
  have hPdeg : P.natDegree = m := by
    rw [hP, Polynomial.natDegree_sub_eq_left_of_natDegree_lt hdeg, Polynomial.natDegree_X_pow]
  have hroot : ∀ s : F, s ^ m = s → s ∈ P.roots := by
    intro s hs
    rw [Polynomial.mem_roots hPmonic.ne_zero, Polynomial.IsRoot.def, hP]
    simp [hs]
  have hsub : T ⊆ P.roots.toFinset := by
    intro s hs
    rw [Multiset.mem_toFinset]
    exact hroot s (hT s hs)
  have hcard : P.roots.toFinset.card ≤ m :=
    le_trans (Multiset.toFinset_card_le _) (le_trans (Polynomial.card_roots' P) (le_of_eq hPdeg))
  have heq : T = P.roots.toFinset := Finset.eq_of_subset_of_card_le hsub (by omega)
  rw [heq, Multiset.mem_toFinset]
  exact hroot t ht

/-- **★★★★★★★★★★★★★★★★`𝓀_{K^ur} → 𝓀_{K̄}` は全射**——`t ∈ 𝓀_{K̄}` は
`t^{q^n} = t` を満たし、`𝓀_{K^ur}` の中の `q^n` 元の部分体がその方程式の根を
**すべて**尽くすから(根はたかだか `q^n` 個)。

★これが「`K̄/K^ur` は完全分岐」の内容だが、惰性群を一切経由していない。 -/
theorem unramResidueToClosure_surjective (K : PAdicLocalField p) :
    Function.Surjective (unramResidueToClosure K) := by
  classical
  intro t
  obtain ⟨n, hn, htn⟩ := exists_pow_eq_self_closureResidueField K t
  set x : K.closure := unramLevelValGen K n with hx
  have hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure) ≤ unramifiedClosure K :=
    adjoin_le_unramifiedClosure K (isUnramified_unramLevelValGen K hn)
  haveI : Fintype (IsLocalRing.ResidueField (adjoinIntegers K x)) := Fintype.ofFinite _
  have hcard : Fintype.card (IsLocalRing.ResidueField (adjoinIntegers K x))
      = Nat.card 𝓀[K.carrier] ^ n := by
    rw [← Nat.card_eq_fintype_card, hx]; exact card_residueField_unramLevel K hn
  set φ := adjResidueToClosure K x with hφ
  set T : Finset (closureResidueField K) := Finset.image φ Finset.univ with hT
  have hTcard : T.card = Nat.card 𝓀[K.carrier] ^ n := by
    rw [hT, Finset.card_image_of_injective _ φ.injective, Finset.card_univ, hcard]
  have hTpow : ∀ s ∈ T, s ^ (Nat.card 𝓀[K.carrier] ^ n) = s := by
    intro s hs
    rw [hT, Finset.mem_image] at hs
    obtain ⟨u, -, rfl⟩ := hs
    rw [← map_pow, ← hcard, FiniteField.pow_card]
  have hmem : t ∈ T :=
    mem_of_pow_eq_self_of_card_eq (Nat.one_lt_pow hn (one_lt_residueCard K)) hTcard hTpow htn
  rw [hT, Finset.mem_image] at hmem
  obtain ⟨u, -, rfl⟩ := hmem
  exact ⟨adjResidueToUnram K hle u, unramResidueToClosure_adjResidueToUnram K hle u⟩

/-- **`𝓀_{K^ur} ≃+* 𝓀_{K̄}`**。単射は体からの環準同型だから、全射は上の数え上げ。 -/
noncomputable def unramResidueFieldEquivClosure (K : PAdicLocalField p) :
    unramResidueField K ≃+* closureResidueField K :=
  RingEquiv.ofBijective (unramResidueToClosure K)
    ⟨(unramResidueToClosure K).injective, unramResidueToClosure_surjective K⟩

/-- **★★★★★★★★★★★★★★★★★★★★(目標 2)`𝓀_{K^ur}` は代数閉体**——すなわち `𝔽̄_q`。

退化の自己検査:`K^ur` を有限次不分岐拡大 `K_N` に替えると**偽**
(`𝓀_{K_N} = 𝔽_{q^N}` は代数閉でない)。無限に登る塔が要る。 -/
theorem isAlgClosed_unramResidueField (K : PAdicLocalField p) :
    IsAlgClosed (unramResidueField K) := by
  haveI := isAlgClosed_residueField_closureInt K
  exact IsAlgClosed.of_ringEquiv _ _ (unramResidueFieldEquivClosure K).symm

/-- **`𝓀_{K^ur}` では `n` 乗が全射**(Dwork の補題の第 0 段)。 -/
theorem exists_pow_eq (K : PAdicLocalField p) (a : unramResidueField K) {n : ℕ} (hn : 0 < n) :
    ∃ z : unramResidueField K, z ^ n = a := by
  haveI := isAlgClosed_unramResidueField K
  exact IsAlgClosed.exists_pow_nat_eq a hn

/-- **Artin–Schreier(代数閉体版)**——`c ↦ c^q - c` は全射(`q > 1`)。

退化の自己検査:`q = 1` では `c^1 - c = 0` なので `a ≠ 0` に対し**偽**。 -/
theorem exists_pow_sub_self_eq_of_isAlgClosed {F : Type*} [Field F] [IsAlgClosed F] {q : ℕ}
    (hq1 : 1 < q) (a : F) : ∃ c : F, c ^ q - c = a := by
  have hXa : (Polynomial.X + Polynomial.C a : Polynomial F).natDegree ≤ 1 := by
    refine le_trans (Polynomial.natDegree_add_le _ _) ?_
    simp
  have hXad : (Polynomial.X + Polynomial.C a : Polynomial F).degree < (q : ℕ) := by
    refine lt_of_le_of_lt Polynomial.degree_le_natDegree ?_
    exact_mod_cast lt_of_le_of_lt hXa hq1
  set P : Polynomial F := Polynomial.X ^ q - (Polynomial.X + Polynomial.C a) with hP
  have hPmonic : P.Monic := Polynomial.monic_X_pow_sub hXad
  have hPnat : P.natDegree = q := by
    rw [hP, Polynomial.natDegree_sub_eq_left_of_natDegree_lt
      (by rw [Polynomial.natDegree_X_pow]; omega), Polynomial.natDegree_X_pow]
  have hPdeg : P.degree ≠ 0 := by
    rw [Polynomial.degree_eq_natDegree hPmonic.ne_zero, hPnat]
    exact_mod_cast (by omega : q ≠ 0)
  obtain ⟨c, hc⟩ := IsAlgClosed.exists_root P hPdeg
  refine ⟨c, ?_⟩
  rw [Polynomial.IsRoot.def, hP] at hc
  simp only [Polynomial.eval_sub, Polynomial.eval_pow, Polynomial.eval_X, Polynomial.eval_add,
    Polynomial.eval_C] at hc
  linear_combination hc

/-- **`𝓀_{K^ur}` では `c ↦ c^q - c` が全射**(Dwork の補題の逐次近似の各段)。 -/
theorem exists_pow_card_sub_self_eq (K : PAdicLocalField p) (a : unramResidueField K) :
    ∃ c : unramResidueField K, c ^ (Nat.card 𝓀[K.carrier]) - c = a := by
  haveI := isAlgClosed_unramResidueField K
  exact exists_pow_sub_self_eq_of_isAlgClosed (one_lt_residueCard K) a

/-! ## 8. 目標 3 —— 算術 Frobenius

★`coherentFrobenius`(`UnramifiedZhat.lean`)は位相的生成元一般であって
算術 Frobenius とは限らない(モジュール docstring の訂正)。ここでは
**剰余体の条件を直接課した**集合のコンパクト性で算術 Frobenius を作る。 -/

/-- `Γ_K` の元は `K.closure` のスペクトルノルムを保つ。 -/
theorem norm_absGal (K : PAdicLocalField p) (g : K.closure ≃ₐ[K.carrier] K.closure)
    (z : K.closure) : ‖g z‖ = ‖z‖ := by
  rw [NormedAlgebra.norm_eq_spectralNorm K.carrier (g z),
    NormedAlgebra.norm_eq_spectralNorm K.carrier z]
  exact (spectralNorm_eq_of_equiv g z).symm

/-- `Γ_K` の元の `𝒪_{K̄}` への制限。 -/
noncomputable def absGalInt (K : PAdicLocalField p) (g : K.closure ≃ₐ[K.carrier] K.closure) :
    ↥(closureInt K) ≃+* ↥(closureInt K) where
  toFun z := ⟨g (z : K.closure), by rw [mem_closureInt, norm_absGal]; exact z.2⟩
  invFun z := ⟨g.symm (z : K.closure), by rw [mem_closureInt, norm_absGal]; exact z.2⟩
  left_inv z := Subtype.ext (g.symm_apply_apply _)
  right_inv z := Subtype.ext (g.apply_symm_apply _)
  map_mul' a b := Subtype.ext (map_mul g _ _)
  map_add' a b := Subtype.ext (map_add g _ _)

@[simp] theorem absGalInt_coe (K : PAdicLocalField p) (g : K.closure ≃ₐ[K.carrier] K.closure)
    (z : ↥(closureInt K)) :
    ((absGalInt K g z : ↥(closureInt K)) : K.closure) = g (z : K.closure) := rfl

/-- `Γ_K` の元が誘導する `𝓀_{K̄}` の自己同型。 -/
noncomputable def absGalResidue (K : PAdicLocalField p)
    (g : K.closure ≃ₐ[K.carrier] K.closure) :
    closureResidueField K ≃+* closureResidueField K :=
  IsLocalRing.ResidueField.mapEquiv (absGalInt K g)

theorem absGalResidue_residue (K : PAdicLocalField p) (g : K.closure ≃ₐ[K.carrier] K.closure)
    (z : ↥(closureInt K)) :
    absGalResidue K g (IsLocalRing.residue _ z)
      = IsLocalRing.residue _ (absGalInt K g z) := by
  rw [absGalResidue, IsLocalRing.ResidueField.mapEquiv_apply, IsLocalRing.ResidueField.map_residue]
  rfl

/-- 開部分群による右移動で閉じた集合は開。 -/
theorem isOpen_of_mul_mem {G : Type*} [Group G] [TopologicalSpace G] [IsTopologicalGroup G]
    {H : Subgroup G} (hH : IsOpen (H : Set G)) {S : Set G}
    (hS : ∀ g ∈ S, ∀ h ∈ H, g * h ∈ S) : IsOpen S := by
  rw [isOpen_iff_forall_mem_open]
  intro g hg
  refine ⟨(fun y => g * y) '' (H : Set G), ?_, (isOpenMap_mul_left g) _ hH,
    ⟨1, one_mem _, mul_one g⟩⟩
  rintro _ ⟨h, hh, rfl⟩
  exact hS g hg h hh

/-- 開部分群による右移動で閉じた集合は閉(補集合も同じ性質を持つから)。 -/
theorem isClosed_of_mul_mem {G : Type*} [Group G] [TopologicalSpace G] [IsTopologicalGroup G]
    {H : Subgroup G} (hH : IsOpen (H : Set G)) {S : Set G}
    (hS : ∀ g ∈ S, ∀ h ∈ H, g * h ∈ S) : IsClosed S := by
  rw [← isOpen_compl_iff]
  refine isOpen_of_mul_mem hH ?_
  intro g hg h hh hcon
  refine hg ?_
  have := hS (g * h) hcon h⁻¹ (inv_mem hh)
  rwa [mul_assoc, mul_inv_cancel, mul_one] at this

/-- `K(x)` の整数の剰余に `q` 乗として作用する `Γ_K` の元の全体。 -/
def absGalFrobSet (K : PAdicLocalField p) (x : K.closure) :
    Set (K.closure ≃ₐ[K.carrier] K.closure) :=
  {g | ∀ z : ↥(closureInt K),
    (z : K.closure) ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure) →
      absGalResidue K g (IsLocalRing.residue _ z)
        = (IsLocalRing.residue _ z) ^ (Nat.card 𝓀[K.carrier])}

theorem absGalFrobSet_mul_mem (K : PAdicLocalField p) (x : K.closure)
    {g : K.closure ≃ₐ[K.carrier] K.closure} (hg : g ∈ absGalFrobSet K x)
    {h : K.closure ≃ₐ[K.carrier] K.closure}
    (hh : h ∈ (IntermediateField.adjoin K.carrier ({x} : Set K.closure)).fixingSubgroup) :
    g * h ∈ absGalFrobSet K x := by
  intro z hz
  have hfix : h (z : K.closure) = (z : K.closure) :=
    (IntermediateField.mem_fixingSubgroup_iff _ h).mp hh _ hz
  have heq : absGalInt K (g * h) z = absGalInt K g z := by
    apply Subtype.ext
    rw [absGalInt_coe, absGalInt_coe, AlgEquiv.mul_apply, hfix]
  rw [absGalResidue_residue, heq, ← absGalResidue_residue]
  exact hg z hz

theorem isClosed_absGalFrobSet (K : PAdicLocalField p) (x : K.closure) :
    IsClosed (absGalFrobSet K x) :=
  isClosed_of_mul_mem
    (IntermediateField.fixingSubgroup_isOpen
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
    (fun _ hg _ hh => absGalFrobSet_mul_mem K x hg hh)

theorem absGalFrobSet_mono (K : PAdicLocalField p) {x y : K.closure}
    (hxy : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({y} : Set K.closure)) :
    absGalFrobSet K y ⊆ absGalFrobSet K x :=
  fun _ hg z hz => hg z (hxy hz)

/-- **`K(x)` の Frobenius は `Γ_K` へ持ち上がる**——`K(x)/K` は normal なので
制限射 `Γ_K ↠ Gal(K(x)/K)` は全射(`AlgEquiv.restrictNormalHom_surjective`)。 -/
theorem absGalFrobSet_nonempty_of_frobenius (K : PAdicLocalField p) (x : K.closure)
    (hu : IsUnramifiedAdjoin K x)
    (σ : ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
    (hσ : ∀ z : IsLocalRing.ResidueField (adjoinIntegers K x),
      residueAlgEquiv K x σ z = z ^ (Nat.card 𝓀[K.carrier])) :
    (absGalFrobSet K x).Nonempty := by
  haveI : Normal K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
    normal_of_isUnramifiedAdjoin K x hu
  obtain ⟨g, hg⟩ := AlgEquiv.restrictNormalHom_surjective
    (F := K.carrier) (K₁ := ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
    K.closure σ
  have hg' : g.restrictNormal
      ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure)) = σ := hg
  refine ⟨g, ?_⟩
  intro z hz
  set z' : ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
    ⟨(z : K.closure), hz⟩ with hz'
  have hnorm : z' ∈ adjoinIntegers K x := (mem_closureInt K _).mp z.2
  set z'' : adjoinIntegers K x := ⟨z', hnorm⟩ with hz''
  have hgz : g (z : K.closure) = ((σ z' : ↥(IntermediateField.adjoin K.carrier
      ({x} : Set K.closure))) : K.closure) := by
    have hcom := AlgEquiv.restrictNormal_commutes g
      (↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure))) z'
    rw [hg'] at hcom
    exact hcom.symm
  have hval : ((algEquivIntegers K x σ z'' : adjoinIntegers K x) :
      ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure))) = σ z' := rfl
  have heq : absGalInt K g z = adjToClosureInt K x (algEquivIntegers K x σ z'') := by
    apply Subtype.ext
    rw [absGalInt_coe, adjToClosureInt_coe, hval]
    exact hgz
  have hzz : adjToClosureInt K x z'' = z := by
    apply Subtype.ext
    rw [adjToClosureInt_coe]
  calc absGalResidue K g (IsLocalRing.residue _ z)
      = IsLocalRing.residue _ (absGalInt K g z) := absGalResidue_residue K g z
    _ = IsLocalRing.residue _ (adjToClosureInt K x (algEquivIntegers K x σ z'')) := by rw [heq]
    _ = adjResidueToClosure K x (IsLocalRing.residue _ (algEquivIntegers K x σ z'')) :=
        (adjResidueToClosure_residue K x _).symm
    _ = adjResidueToClosure K x (residueAlgEquiv K x σ (IsLocalRing.residue _ z'')) := by
        rw [residueAlgEquiv_apply]
    _ = adjResidueToClosure K x ((IsLocalRing.residue _ z'') ^ (Nat.card 𝓀[K.carrier])) := by
        rw [hσ]
    _ = (adjResidueToClosure K x (IsLocalRing.residue _ z'')) ^ (Nat.card 𝓀[K.carrier]) :=
        map_pow _ _ _
    _ = (IsLocalRing.residue _ (adjToClosureInt K x z'')) ^ (Nat.card 𝓀[K.carrier]) := by
        rw [adjResidueToClosure_residue]
    _ = (IsLocalRing.residue _ z) ^ (Nat.card 𝓀[K.carrier]) := by rw [hzz]

/-- 次数 `N` の不分岐拡大の「Frobenius つき生成元」。`N = 0` では `0`(意味を持たない)。 -/
noncomputable def frobGen (K : PAdicLocalField p) (N : ℕ) : K.closure :=
  if h : N = 0 then 0 else (exists_frobenius K N h).choose

theorem frobGen_spec (K : PAdicLocalField p) {N : ℕ} (hN : N ≠ 0) :
    Module.finrank K.carrier
        (IntermediateField.adjoin K.carrier ({frobGen K N} : Set K.closure)) = N
      ∧ IsUnramifiedAdjoin K (frobGen K N)
      ∧ ∃ σ : ↥(IntermediateField.adjoin K.carrier ({frobGen K N} : Set K.closure))
          ≃ₐ[K.carrier] ↥(IntermediateField.adjoin K.carrier ({frobGen K N} : Set K.closure)),
        (∀ z : IsLocalRing.ResidueField (adjoinIntegers K (frobGen K N)),
            residueAlgEquiv K (frobGen K N) σ z = z ^ (Nat.card 𝓀[K.carrier]))
          ∧ orderOf σ = N := by
  rw [frobGen, dif_neg hN]
  exact (exists_frobenius K N hN).choose_spec

theorem absGalFrobSet_frobGen_nonempty (K : PAdicLocalField p) {N : ℕ} (hN : N ≠ 0) :
    (absGalFrobSet K (frobGen K N)).Nonempty := by
  obtain ⟨-, hu, σ, hσ, -⟩ := frobGen_spec K hN
  exact absGalFrobSet_nonempty_of_frobenius K _ hu σ hσ

theorem absGalFrobSet_frobGen_mono (K : PAdicLocalField p) {M N : ℕ} (hM : M ≠ 0) (hN : N ≠ 0)
    (h : M ∣ N) : absGalFrobSet K (frobGen K N) ⊆ absGalFrobSet K (frobGen K M) :=
  absGalFrobSet_mono K
    (adjoin_le_of_dvd K _ _ (frobGen_spec K hM).2.1 (frobGen_spec K hN).2.1
      (by rw [(frobGen_spec K hM).1, (frobGen_spec K hN).1]; exact h))

/-- **すべての段で Frobenius として作用する `Γ_K` の元が存在する**(コンパクト性)。

各段の集合は空でない閉集合(`Γ_K` はコンパクト)で、`N ∣ NM` により有向に減少する。 -/
theorem exists_absGal_frobenius (K : PAdicLocalField p) :
    ∃ g : K.closure ≃ₐ[K.carrier] K.closure,
      ∀ N : ℕ, N ≠ 0 → g ∈ absGalFrobSet K (frobGen K N) := by
  haveI : Nonempty {N : ℕ // N ≠ 0} := ⟨⟨1, one_ne_zero⟩⟩
  obtain ⟨g, hg⟩ := IsCompact.nonempty_iInter_of_directed_nonempty_isCompact_isClosed
    (fun i : {N : ℕ // N ≠ 0} => absGalFrobSet K (frobGen K i.1))
    (fun i j => ⟨⟨i.1 * j.1, mul_ne_zero i.2 j.2⟩,
      absGalFrobSet_frobGen_mono K i.2 (mul_ne_zero i.2 j.2) ⟨j.1, rfl⟩,
      absGalFrobSet_frobGen_mono K j.2 (mul_ne_zero i.2 j.2) ⟨i.1, mul_comm i.1 j.1⟩⟩)
    (fun i => absGalFrobSet_frobGen_nonempty K i.2)
    (fun i => (isClosed_absGalFrobSet K _).isCompact)
    (fun i => isClosed_absGalFrobSet K _)
  exact ⟨g, fun N hN => Set.mem_iInter.mp hg ⟨N, hN⟩⟩

/-- **`Γ_K` の Frobenius**——`𝒪_{K̄}` の `K^ur` に入る元の剰余に `q` 乗として作用する。

`K^ur` の元はどれかの不分岐 `K(y)` に入り(`mem_unramifiedClosure_iff`)、
その `K(y)` は同じ次数の `K(frobGen K n)` に入る(`adjoin_le_of_dvd`)。 -/
theorem exists_absGal_frobenius_unramified (K : PAdicLocalField p) :
    ∃ g : K.closure ≃ₐ[K.carrier] K.closure, ∀ z : ↥(closureInt K),
      (z : K.closure) ∈ unramifiedClosure K →
        absGalResidue K g (IsLocalRing.residue _ z)
          = (IsLocalRing.residue _ z) ^ (Nat.card 𝓀[K.carrier]) := by
  obtain ⟨g, hg⟩ := exists_absGal_frobenius K
  refine ⟨g, fun z hz => ?_⟩
  obtain ⟨y, hy, hzy⟩ := (mem_unramifiedClosure_iff K (z : K.closure)).mp hz
  have hn : Module.finrank K.carrier
      (IntermediateField.adjoin K.carrier ({y} : Set K.closure)) ≠ 0 := Module.finrank_pos.ne'
  have hle : IntermediateField.adjoin K.carrier ({y} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({frobGen K (Module.finrank K.carrier
          (IntermediateField.adjoin K.carrier ({y} : Set K.closure)))} : Set K.closure) :=
    adjoin_le_of_dvd K y _ hy (frobGen_spec K hn).2.1 (by rw [(frobGen_spec K hn).1])
  exact hg _ hn z (hle hzy)

/-- `Gal(K^ur/K)` の元の `𝒪_{K^ur}` への制限。 -/
noncomputable def unramGalInt (K : PAdicLocalField p) (σ : unramGal K) :
    ↥(unramifiedClosureInt K) ≃+* ↥(unramifiedClosureInt K) where
  toFun w := ⟨σ (w : ↥(unramifiedClosure K)), by
    rw [mem_unramifiedClosureInt, norm_unramGal]
    exact (mem_unramifiedClosureInt K _).mp w.2⟩
  invFun w := ⟨σ.symm (w : ↥(unramifiedClosure K)), by
    rw [mem_unramifiedClosureInt, norm_unramGal]
    exact (mem_unramifiedClosureInt K _).mp w.2⟩
  left_inv w := Subtype.ext (σ.symm_apply_apply _)
  right_inv w := Subtype.ext (σ.apply_symm_apply _)
  map_mul' a b := Subtype.ext (map_mul σ _ _)
  map_add' a b := Subtype.ext (map_add σ _ _)

theorem unramGalInt_apply (K : PAdicLocalField p) (σ : unramGal K)
    (w : ↥(unramifiedClosureInt K)) :
    unramGalInt K σ w = ⟨σ (w : ↥(unramifiedClosure K)), by
      rw [mem_unramifiedClosureInt, norm_unramGal]
      exact (mem_unramifiedClosureInt K _).mp w.2⟩ := rfl

theorem unramGalInt_coe (K : PAdicLocalField p) (σ : unramGal K)
    (w : ↥(unramifiedClosureInt K)) :
    ((unramGalInt K σ w : ↥(unramifiedClosureInt K)) : ↥(unramifiedClosure K))
      = σ (w : ↥(unramifiedClosure K)) := by
  rw [unramGalInt_apply]

/-- `Gal(K^ur/K)` の元が誘導する `𝓀_{K^ur}` の自己同型。 -/
noncomputable def unramGalResidue (K : PAdicLocalField p) (σ : unramGal K) :
    unramResidueField K ≃+* unramResidueField K :=
  IsLocalRing.ResidueField.mapEquiv (unramGalInt K σ)

theorem unramGalResidue_residue (K : PAdicLocalField p) (σ : unramGal K)
    (w : ↥(unramifiedClosureInt K)) :
    unramGalResidue K σ (IsLocalRing.residue _ w)
      = IsLocalRing.residue _ (unramGalInt K σ w) := by
  rw [unramGalResidue, IsLocalRing.ResidueField.mapEquiv_apply,
    IsLocalRing.ResidueField.map_residue]
  rfl

/-- **★★★★★★★★★★★★★★★★★★★★★★(目標 3)算術 Frobenius の存在**——
`Gal(K^ur/K)` には、`𝓀_{K^ur}` 上で `z ↦ z^q` として作用する元がある。

★`UnramifiedZhat.lean` の `coherentFrobenius` は**位相的生成元にすぎず**、
剰余体上の作用は一般に `z ↦ z^{q^k}`(`k ∈ Ẑ^×`)なので**これとは違う**
(モジュール docstring の訂正)。

退化の自己検査。

* `q` を他の数に替えると**偽**——`Gal(𝓀_{K^ur}/𝓀_K)` の中で
  `z ↦ z^q` は各段で位数がちょうど段の次数になる特別な元。
* `σ` を `∃` の外に出すと自由なパラメータになるので内側に閉じ込めてある。 -/
theorem exists_arithFrobenius (K : PAdicLocalField p) :
    ∃ σ : unramGal K, ∀ w : ↥(unramifiedClosureInt K),
      unramGalResidue K σ (IsLocalRing.residue _ w)
        = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier]) := by
  haveI : Normal K.carrier ↥(unramifiedClosure K) := normal_unramifiedClosure K
  obtain ⟨g, hg⟩ := exists_absGal_frobenius_unramified K
  refine ⟨AlgEquiv.restrictNormalHom ↥(unramifiedClosure K) g, fun w => ?_⟩
  have hcom : (((AlgEquiv.restrictNormalHom (F := K.carrier) ↥(unramifiedClosure K) g)
        (w : ↥(unramifiedClosure K)) : ↥(unramifiedClosure K)) : K.closure)
      = g ((w : ↥(unramifiedClosure K)) : K.closure) :=
    AlgEquiv.restrictNormal_commutes g (↥(unramifiedClosure K)) (w : ↥(unramifiedClosure K))
  have hcompat : urToClosureInt K
        (unramGalInt K (AlgEquiv.restrictNormalHom ↥(unramifiedClosure K) g) w)
      = absGalInt K g (urToClosureInt K w) := by
    apply Subtype.ext
    rw [urToClosureInt_coe, absGalInt_coe, urToClosureInt_coe, unramGalInt_coe]
    exact hcom
  have hmemur : ((urToClosureInt K w : ↥(closureInt K)) : K.closure) ∈ unramifiedClosure K := by
    rw [urToClosureInt_coe]
    exact (w : ↥(unramifiedClosure K)).2
  refine (unramResidueToClosure K).injective ?_
  rw [unramGalResidue_residue, unramResidueToClosure_residue, hcompat,
    ← absGalResidue_residue, ← unramResidueToClosure_residue, map_pow]
  exact hg _ hmemur

/-! ## 9. Λ6(Dwork の補題)の入口

完備化側へ渡す。`UnramifiedCompletion.residueFieldEquiv`(完備化は剰余体を
変えない)で全部運べる。 -/

/-- `𝓀_{K̂^ur}` も代数閉体——完備化は剰余体を変えない(`residueFieldEquiv`)。 -/
theorem isAlgClosed_residueField_unramifiedCompletionInt (K : PAdicLocalField p) :
    IsAlgClosed (IsLocalRing.ResidueField ↥(unramifiedCompletionInt K)) := by
  haveI := isAlgClosed_unramResidueField K
  exact IsAlgClosed.of_ringEquiv _ _ (residueFieldEquiv K)

/-- `𝒪_{K^ur} → 𝒪_{K̂^ur}` は Galois 作用と可換。 -/
theorem unramifiedIntHom_unramGalInt (K : PAdicLocalField p) (σ : unramGal K)
    (w : ↥(unramifiedClosureInt K)) :
    unramifiedIntHom K (unramGalInt K σ w)
      = unramGalCompletionInt K σ (unramifiedIntHom K w) := by
  apply Subtype.ext
  rw [unramifiedIntHom_coe, unramGalCompletionInt_coe, unramifiedIntHom_coe, unramGalInt_coe,
    unramGalCompletion_coe]

/-- `Gal(K^ur/K)` の元が誘導する `𝓀_{K̂^ur}` の自己同型。 -/
noncomputable def unramGalCompletionResidue (K : PAdicLocalField p) (σ : unramGal K) :
    IsLocalRing.ResidueField ↥(unramifiedCompletionInt K)
      ≃+* IsLocalRing.ResidueField ↥(unramifiedCompletionInt K) :=
  IsLocalRing.ResidueField.mapEquiv (unramGalCompletionInt K σ)

theorem unramGalCompletionResidue_residue (K : PAdicLocalField p) (σ : unramGal K)
    (v : ↥(unramifiedCompletionInt K)) :
    unramGalCompletionResidue K σ (IsLocalRing.residue _ v)
      = IsLocalRing.residue _ (unramGalCompletionInt K σ v) := by
  rw [unramGalCompletionResidue, IsLocalRing.ResidueField.mapEquiv_apply,
    IsLocalRing.ResidueField.map_residue]
  rfl

theorem residueFieldEquiv_residue (K : PAdicLocalField p) (w : ↥(unramifiedClosureInt K)) :
    residueFieldEquiv K (IsLocalRing.residue _ w)
      = IsLocalRing.residue _ (unramifiedIntHom K w) :=
  IsLocalRing.ResidueField.map_residue _ w

theorem residueFieldEquiv_unramGalResidue (K : PAdicLocalField p) (σ : unramGal K)
    (t : unramResidueField K) :
    residueFieldEquiv K (unramGalResidue K σ t)
      = unramGalCompletionResidue K σ (residueFieldEquiv K t) := by
  obtain ⟨w, rfl⟩ := IsLocalRing.residue_surjective t
  rw [unramGalResidue_residue, residueFieldEquiv_residue, residueFieldEquiv_residue,
    unramGalCompletionResidue_residue, unramifiedIntHom_unramGalInt]

/-- **★★★★★★★★★★★★★★★★★★★★★★★★(Λ6 の入口)`K̂^ur` の算術 Frobenius**——
`Gal(K^ur/K)` には、その `𝒪_{K̂^ur}` への延長が剰余体上で `z ↦ z^q` として
作用する元がある。

これと `exists_pow_eq_completion` / `exists_pow_card_sub_self_eq_completion`、
`UnramifiedCompletion.exists_norm_le_one_norm_sub_lt`(整数の任意精度近似)が
Dwork の補題(`Φ(ξ)/ξ = u`)の逐次近似の材料一式である。 -/
theorem exists_arithFrobenius_completion (K : PAdicLocalField p) :
    ∃ σ : unramGal K, ∀ v : ↥(unramifiedCompletionInt K),
      IsLocalRing.residue _ (unramGalCompletionInt K σ v)
        = (IsLocalRing.residue _ v) ^ (Nat.card 𝓀[K.carrier]) := by
  obtain ⟨σ, hσ⟩ := exists_arithFrobenius K
  refine ⟨σ, fun v => ?_⟩
  obtain ⟨t, ht⟩ := (residueFieldEquiv K).surjective (IsLocalRing.residue _ v)
  calc IsLocalRing.residue _ (unramGalCompletionInt K σ v)
      = unramGalCompletionResidue K σ (IsLocalRing.residue _ v) :=
        (unramGalCompletionResidue_residue K σ v).symm
    _ = unramGalCompletionResidue K σ (residueFieldEquiv K t) := by rw [ht]
    _ = residueFieldEquiv K (unramGalResidue K σ t) :=
        (residueFieldEquiv_unramGalResidue K σ t).symm
    _ = residueFieldEquiv K (t ^ (Nat.card 𝓀[K.carrier])) := by
        obtain ⟨w, rfl⟩ := IsLocalRing.residue_surjective t
        rw [hσ w]
    _ = (residueFieldEquiv K t) ^ (Nat.card 𝓀[K.carrier]) := map_pow _ _ _
    _ = (IsLocalRing.residue _ v) ^ (Nat.card 𝓀[K.carrier]) := by rw [ht]

/-- `𝓀_{K̂^ur}` で `n` 乗は全射。 -/
theorem exists_pow_eq_completion (K : PAdicLocalField p)
    (a : IsLocalRing.ResidueField ↥(unramifiedCompletionInt K)) {n : ℕ} (hn : 0 < n) :
    ∃ z : IsLocalRing.ResidueField ↥(unramifiedCompletionInt K), z ^ n = a := by
  haveI := isAlgClosed_residueField_unramifiedCompletionInt K
  exact IsAlgClosed.exists_pow_nat_eq a hn

/-- `𝓀_{K̂^ur}` で `c ↦ c^q - c` は全射(Dwork の補題の逐次近似の各段)。 -/
theorem exists_pow_card_sub_self_eq_completion (K : PAdicLocalField p)
    (a : IsLocalRing.ResidueField ↥(unramifiedCompletionInt K)) :
    ∃ c : IsLocalRing.ResidueField ↥(unramifiedCompletionInt K),
      c ^ (Nat.card 𝓀[K.carrier]) - c = a := by
  haveI := isAlgClosed_residueField_unramifiedCompletionInt K
  exact exists_pow_sub_self_eq_of_isAlgClosed (one_lt_residueCard K) a

end ABC3.Found.PGC
