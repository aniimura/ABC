import ABC3.Found.FrdI.Def31Pf

/-!
# [FrdI] Proposition 3.2 —— Frobenioid の perfection

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.58–p.59。

原文 (FrdI p.58):
> (Perfections of Frobenioids) Suppose that the Frobenioid

## ★主張の規模(測定)

`Proposition 3.2` は **3 条**:

| 条 | 内容 | 状態 |
|---|---|---|
| (i) | 1-可換図式 `𝒞 → 𝒞^pf` / `𝒞 → 𝔽_Φ` / `𝒞^pf → 𝔽_{Φ^pf}` / `𝔽_Φ → 𝔽_{Φ^pf}`、および `𝒞^pf` の pre-Frobenioid 構造 | ★**ここで実装** |
| (ii) | `𝒞^pf` の射の 10 種類の判定を「cofinal collection」で述べる | 未 |
| (iii) | `𝒞^pf` が perfect かつ isotropic 型の Frobenioid、`𝒞^pf ≃ (𝒞^pf)^pf` | 未 |

## ★★(i) の急所 —— 「次数で割る」

`f ∈ Hom^pf(A,B)` を `(a : A→A′, b : B→B′)`(次数 `n`)と `φ : A′→B′` で表すとき、
`𝔽_{Φ^pf}` へ落とす 3 つの量は:

| 量 | 値 | 代表元の取り替えで不変な理由 |
|---|---|---|
| **次数** | `degFr φ` | `Proposition 1.10, (i)` の `degFr(φ) = degFr(φ′)` |
| **底** | `Base a ≫ Base φ ≫ (Base b)⁻¹` | `Base` の関手性と `a`・`b` の base-iso 性 |
| **零因子** | `(Base a)^*(Div φ) / n ∈ Φ^pf(A_𝒟)` | ★`Div(φ′) = degFr(β)·(Base α)⁻¹^*(Div φ)` と `Pf` の同一視 |

★★**零因子だけが `Φ` では書けず `Φ^pf` が要る** —— これが「perfection」の名の由来である。

## ★★宇宙について(実装上の判断)

★`Hom^pf` は `Type (max u2 v2)` に住むが、`𝒟` の射は `Type v`、`Φ^pf(A_𝒟)` は `Type w` に住む。
★**余錐(`colimit.desc`)で降ろすと宇宙が合わない。**
★★**そこで `HomPf.exists_rep` で代表元を選び(`Classical.choose`)、
「代表元の取り替えで不変」を別に示す**という形にする。
★これなら宇宙の制約が一切かからない(`Def31Pf.lean` を書き換えずに済む)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory CategoryTheory.Limits

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (F : FrobenioidCore P)

/-! ## ★0. `Pf` の道具 —— **分母が共通なら分母は増えない**

★★`Pf.mk_add_mk` は `mk a m + mk b n = mk (n•a + m•b) (m*n)` なので、
**`m = n` のときも分母を掛けてしまう**。★そのまま零因子の合成則に使うと
`(n*n) • (a+b) = n • (n•a + n•b)` という係数の突き合わせが目標の**内側**に残り、
`Φ.map` の coercion に阻まれて `simp` が当たらない(2026-08-17 に測定)。

★★**そこで「分母共通版」を先に 1 本置く。** `Φ` も `P` も現れない
**純粋に `Pf` だけの主張**なので、係数計算が coercion に邪魔されない。 -/

/-- ★★**分母が共通なら分母は増えない**。 -/
theorem Pf.mk_add_mk_same {M : Type w} [AddCommMonoid M] (m m' : M) (a : ℕ+) :
    Pf.mk m a + Pf.mk m' a = Pf.mk (m + m') a := by
  rw [Pf.mk_add_mk]
  refine Pf.sound 1 ?_
  push_cast
  simp only [one_mul, smul_add, smul_smul]

/-! ## ★0-b. **`Pf M` は divisorial を受け継ぐ**

★★段 D(`𝒞^pf` の pre-Frobenioid 構造)が要求する `(Φ^pf).IsDivisorialOn` のため。
`IsDivisorial = (integral ∧ saturated ∧ of-characteristic-type) ∧ sharp` で、
★**sharp は `Pf.isSharp_pf`、of-characteristic-type は
`isOfCharacteristicType_of_isSharp` で既済**。★残る 2 つをここで示す。

★★**これらは `MonoidVocabulary.lean` に属する汎用補題であり、ここは仮置きである。** -/

/-- ★**簡約性は `M^pf` に遺伝する**(左簡約の本体)。

★`mk a p + mk b q = mk a p + mk c r` を `Quotient.exact` で開くと、
両辺の `a` の係数が `k·p·r·q` と `k·p·q·r` で**一致する**ので、`M` の簡約性で消せる。 -/
theorem Pf.add_left_cancel' {M : Type w} [AddCommMonoid M] [IsCancelAdd M]
    (x y z : Pf M) (h : x + y = x + z) : y = z := by
  induction x using Pf.inductionOn with | _ a p =>
  induction y using Pf.inductionOn with | _ b q =>
  induction z using Pf.inductionOn with | _ c r =>
  rw [Pf.mk_add_mk, Pf.mk_add_mk] at h
  obtain ⟨k, e⟩ := Quotient.exact h
  refine Pf.sound (k * p * p) ?_
  push_cast at e ⊢
  simp only [smul_add, smul_smul] at e
  rw [show (k : ℕ) * ((p : ℕ) * (q : ℕ)) * (r : ℕ)
      = (k : ℕ) * ((p : ℕ) * (r : ℕ)) * (q : ℕ) from by ring] at e
  have h2 := add_left_cancel e
  rw [show (k : ℕ) * (p : ℕ) * (p : ℕ) * (r : ℕ)
        = (k : ℕ) * ((p : ℕ) * (r : ℕ)) * (p : ℕ) from by ring,
      show (k : ℕ) * (p : ℕ) * (p : ℕ) * (q : ℕ)
        = (k : ℕ) * ((p : ℕ) * (q : ℕ)) * (p : ℕ) from by ring]
  exact h2

/-- ★**簡約性は `M^pf` に遺伝する**。 -/
theorem Pf.isCancelAdd' {M : Type w} [AddCommMonoid M] [IsCancelAdd M] :
    IsCancelAdd (Pf M) where
  add_left_cancel x y z h := Pf.add_left_cancel' x y z h
  add_right_cancel x y z h :=
    Pf.add_left_cancel' x y z (by rw [add_comm x y, add_comm x z]; exact h)

/-- ★★**`M^pf` は saturated**(`M` が簡約的なら)。

★★**要点: 「perfect + 簡約的 ⟹ saturated」** —— 完全化そのものが saturated 性を生む。
`n • a` が `M^pf` の像に入るとき、その代表を **`M^pf` の perfect 性で `n` 等分**すると、
分母がそのまま割り切れて `a` 自身が像に入る。 -/
theorem Pf.isSaturatedMonoid' {M : Type w} [AddCommMonoid M] [IsCancelAdd M] :
    IsSaturatedMonoid (Pf M) := by
  letI : IsCancelAdd (Pf M) := Pf.isCancelAdd'
  intro a n hn
  induction a using AddLocalization.induction_on with
  | _ y =>
    intro hmem
    obtain ⟨x, hx⟩ := hmem
    rw [AddLocalization.mk_nsmul, toGp, AddLocalization.mk_eq_mk_iff,
      AddLocalization.r_iff_exists] at hx
    obtain ⟨c, hc⟩ := hx
    have h1 : (n • (y.2 : Pf M)) + x = n • y.1 := by
      have := add_left_cancel hc
      simpa [-Pf.zero_def] using this
    obtain ⟨z, hz⟩ := (Pf.isPerfectMonoid_pf (M := M) ⟨n, hn⟩).2 x
    have hz' : n • z = x := hz
    have h2 : n • ((y.2 : Pf M) + z) = n • y.1 := by
      rw [smul_add, hz']
      exact h1
    have h3 : (y.2 : Pf M) + z = y.1 :=
      (Pf.isPerfectMonoid_pf (M := M) ⟨n, hn⟩).1 h2
    refine ⟨z, ?_⟩
    rw [toGp, AddLocalization.mk_eq_mk_iff, AddLocalization.r_iff_exists]
    exact ⟨0, by simpa [-Pf.zero_def] using h3⟩

/-- ★★★**divisorial 性は `M^pf` に遺伝する**。

★4 条件のうち **sharp と of-characteristic-type は既存の補題そのまま**、
残る integral と saturated が上の 2 本である。 -/
theorem Pf.isDivisorial' {M : Type w} [AddCommMonoid M] (h : IsDivisorial M) :
    IsDivisorial (Pf M) := by
  letI := isCancelAdd_of_isIntegralMonoid M h.1.1
  letI : IsCancelAdd (Pf M) := Pf.isCancelAdd'
  have hsharp : IsSharp (Pf M) := Pf.isSharp_pf h.2
  exact ⟨⟨isIntegralMonoid_of_isCancelAdd (Pf M), Pf.isSaturatedMonoid',
    isOfCharacteristicType_of_isSharp (Pf M) hsharp⟩, hsharp⟩

/-! ## ★0-c. `Pf` の「`k` で割る」

★★段 B′(対象が対 `(A,n)` の `𝒞^pf`)で要る —— `(A,n)` は `A` の `n` 乗根なので、
★**零因子も `n` で割られる**。`Pf` では**分母を `n` 倍する**ことである。 -/

/-- ★★**`k` で割る** —— 分母を `k` 倍する。★`Pf` が perfect なので
これは `k • (−)` の逆写像である。 -/
noncomputable def Pf.divBy {M : Type w} [AddCommMonoid M] (k : ℕ+) (x : Pf M) : Pf M :=
  Quotient.liftOn x (fun p => Pf.mk p.1 (k * p.2)) (by
    rintro ⟨m, a⟩ ⟨m', a'⟩ ⟨j, e⟩
    refine Pf.sound j ?_
    have h := congrArg (fun t : M => (k : ℕ) • t) e
    simp only [smul_smul] at h
    push_cast
    rw [show (j : ℕ) * ((k : ℕ) * (a' : ℕ)) = (k : ℕ) * ((j : ℕ) * (a' : ℕ)) from by ring,
      show (j : ℕ) * ((k : ℕ) * (a : ℕ)) = (k : ℕ) * ((j : ℕ) * (a : ℕ)) from by ring]
    exact h)

@[simp] theorem Pf.divBy_mk {M : Type w} [AddCommMonoid M] (k : ℕ+) (m : M) (a : ℕ+) :
    Pf.divBy k (Pf.mk m a) = Pf.mk m (k * a) := rfl

/-- ★**`k` で割ってから `k` 倍すると元に戻る**。 -/
theorem Pf.nsmul_divBy {M : Type w} [AddCommMonoid M] (k : ℕ+) (x : Pf M) :
    ((k : ℕ+) : ℕ) • Pf.divBy k x = x := by
  induction x using Pf.inductionOn with | _ m a =>
  rw [Pf.divBy_mk, Pf.nsmul_mk]
  refine Pf.sound 1 ?_
  push_cast
  simp [smul_smul, mul_comm, mul_assoc, mul_left_comm]

/-! ## ★1. 代表元が定める 3 つの量 -/

variable {P F} in
/-- ★代表元の**次数**。 -/
noncomputable def repDeg {A B : C} (Z : IdxPf P F A B) (φ : Z.right.obj.1 ⟶ Z.right.obj.2) :
    ℕ+ := P.degFr φ

variable {P F} in
/-- ★代表元の**根の次数**(添字対象の脚の Frobenius 次数)。 -/
noncomputable def repRoot {A B : C} (Z : IdxPf P F A B) : ℕ+ := P.degFr Z.hom.hom.1

/-- ★★**底**(`IsIso` を**仮引数**に取る形)。

★★`Prop44.lean` で学んだ技である —— `Z.hom.hom.2` の形で `haveI` すると
インスタンス探索が失敗するので、**素の射と `IsIso` の仮引数に割る**。 -/
noncomputable def repBaseOf {A B A' B' : C} (a : A ⟶ A') (b : B ⟶ B')
    (_hb : IsIso (P.Base b)) (φ : A' ⟶ B') :
    (P.toElem.obj A).base ⟶ (P.toElem.obj B).base :=
  haveI := _hb
  P.Base a ≫ P.Base φ ≫ inv (P.Base b)

variable {P F} in
/-- ★代表元の**底** `Base(a) ≫ Base(φ) ≫ Base(b)⁻¹`。 -/
noncomputable def repBase {A B : C} (Z : IdxPf P F A B) (φ : Z.right.obj.1 ⟶ Z.right.obj.2) :
    (P.toElem.obj A).base ⟶ (P.toElem.obj B).base :=
  repBaseOf P Z.hom.hom.1 Z.hom.hom.2 Z.hom.property.2.1.2 φ

variable {P F} in
/-- ★★代表元の**零因子**(`Φ^pf` の元)—— `(Base a)^*(Div φ) / n`。 -/
noncomputable def repDiv {A B : C} (Z : IdxPf P F A B) (φ : Z.right.obj.1 ⟶ Z.right.obj.2) :
    Pf (Φ.val (P.toElem.obj A).base) :=
  Pf.mk (Φ.map (P.Base Z.hom.hom.1) (P.Div φ)) (repRoot Z)

/-! ## ★2. 遷移で不変であること -/

variable {P F} in
/-- ★★**次数は遷移で不変** —— `Proposition 1.10, (i)` の `degFr(φ) = degFr(φ′)`。 -/
theorem repDeg_map {A B : C} {Z W : IdxPf P F A B} (u : Z ⟶ W)
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) :
    repDeg W (idxTransport P F u φ) = repDeg Z φ :=
  (prop_1_10_i_degFr_phi_eq P u.right.property.2.2 (idxTransport_spec u φ)).symm

/-- ★★**底は遷移で不変**(素の射で述べた形)。 -/
theorem repBaseOf_map {A B A' B' A'' B'' : C}
    (a : A ⟶ A') (b : B ⟶ B') (hb : IsIso (P.Base b))
    (α : A' ⟶ A'') (β : B' ⟶ B'') (_hβ : IsIso (P.Base β))
    (hbβ : IsIso (P.Base (b ≫ β)))
    (φ : A' ⟶ B') (ψ : A'' ⟶ B'') (hsq : φ ≫ β = α ≫ ψ) :
    repBaseOf P (a ≫ α) (b ≫ β) hbβ ψ = repBaseOf P a b hb φ := by
  haveI := hb
  haveI := hbβ
  show P.Base (a ≫ α) ≫ P.Base ψ ≫ inv (P.Base (b ≫ β))
    = P.Base a ≫ P.Base φ ≫ inv (P.Base b)
  have hsq' := congrArg P.Base hsq
  simp only [P.Base_comp] at hsq'
  refine Eq.trans (Category.assoc _ _ _).symm ?_
  refine (IsIso.comp_inv_eq _).mpr ?_
  simp only [P.Base_comp, Category.assoc, IsIso.inv_hom_id_assoc]
  exact congrArg (fun t => P.Base a ≫ t) hsq'.symm

/-- ★`IsIso` は `Prop` なので、射が等しければ値も等しい。 -/
theorem repBaseOf_congr {A B A' B' : C} {a a' : A ⟶ A'} (ha : a = a') {b b' : B ⟶ B'}
    (hb : b = b') (h1 : IsIso (P.Base b)) (h2 : IsIso (P.Base b')) (φ : A' ⟶ B') :
    repBaseOf P a b h1 φ = repBaseOf P a' b' h2 φ := by
  subst ha
  subst hb
  rfl

variable {P F} in
/-- ★★**底は遷移で不変**。 -/
theorem repBase_map {A B : C} {Z W : IdxPf P F A B} (u : Z ⟶ W)
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) :
    repBase W (idxTransport P F u φ) = repBase Z φ := by
  have hw : Z.hom.hom.1 ≫ u.right.hom.1 = W.hom.hom.1 :=
    congrArg (fun t : biFrObj P F A B ⟶ W.right => t.hom.1) (Under.w u)
  have hw2 : Z.hom.hom.2 ≫ u.right.hom.2 = W.hom.hom.2 :=
    congrArg (fun t : biFrObj P F A B ⟶ W.right => t.hom.2) (Under.w u)
  have hiso : IsIso (P.Base (Z.hom.hom.2 ≫ u.right.hom.2)) := hw2 ▸ W.hom.property.2.1.2
  show repBaseOf P W.hom.hom.1 W.hom.hom.2 W.hom.property.2.1.2 (idxTransport P F u φ) = _
  refine Eq.trans (repBaseOf_congr P hw.symm hw2.symm W.hom.property.2.1.2 hiso _) ?_
  exact repBaseOf_map P Z.hom.hom.1 Z.hom.hom.2 Z.hom.property.2.1.2
    u.right.hom.1 u.right.hom.2 u.right.property.2.1.2 hiso φ _ (idxTransport_spec u φ)

variable {P F} in
/-- ★★★**零因子は遷移で不変**(`Φ^pf` の中で)。

★★**これが `Proposition 3.2` の急所**である ——
`Proposition 1.10, (i)` の `Div(φ′) = degFr(β)·(Base α)^*(Div φ)` が
「分子が `m` 倍、分母(根の次数)も `m` 倍」を与え、
`Pf` の同一視でちょうど打ち消し合う。
★★**`Φ` のままでは打ち消せない —— そこが「perfection が要る」ことの内容である。** -/
theorem repDiv_map {A B : C} {Z W : IdxPf P F A B} (u : Z ⟶ W)
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) :
    repDiv W (idxTransport P F u φ) = repDiv Z φ := by
  have hw : Z.hom.hom.1 ≫ u.right.hom.1 = W.hom.hom.1 :=
    congrArg (fun t : biFrObj P F A B ⟶ W.right => t.hom.1) (Under.w u)
  have hdiv : Φ.map (P.Base u.right.hom.1) (P.Div (idxTransport P F u φ))
      = ((P.degFr u.right.hom.2 : ℕ+) : ℕ) • P.Div φ :=
    prop_1_10_i_Div_formula P φ u.right.hom.1 u.right.hom.2 (idxTransport P F u φ)
      u.right.property.1.1.2 u.right.property.2.1.1.2 (idxTransport_spec u φ)
  have hroot : repRoot W = P.degFr u.right.hom.1 * repRoot Z := by
    show P.degFr W.hom.hom.1 = P.degFr u.right.hom.1 * P.degFr Z.hom.hom.1
    rw [← hw, P.degFr_comp]
  have hval : Φ.map (P.Base W.hom.hom.1) (P.Div (idxTransport P F u φ))
      = ((P.degFr u.right.hom.1 : ℕ+) : ℕ) • Φ.map (P.Base Z.hom.hom.1) (P.Div φ) := by
    rw [← hw, P.Base_comp, Φ.map_comp, hdiv, map_nsmul, u.right.property.2.2]
  show Pf.mk (Φ.map (P.Base W.hom.hom.1) (P.Div (idxTransport P F u φ))) (repRoot W)
    = Pf.mk (Φ.map (P.Base Z.hom.hom.1) (P.Div φ)) (repRoot Z)
  rw [hval, hroot]
  refine Pf.sound 1 ?_
  push_cast
  simp only [one_mul, smul_smul]
  exact congrArg (fun t : ℕ => t • Φ.map (P.Base Z.hom.hom.1) (P.Div φ)) (mul_comm _ _)

/-! ## ★3. 代表元によらないこと

★★添字圏が **filtered** なので、2 つの代表元は共通の上界へ持ち上がる。
そこで一致するので、上の 3 つの不変性から**値が代表元によらない**。 -/

variable {P F} in
theorem repDeg_congr {A B : C} {Z Z' : IdxPf P F A B}
    {φ : Z.right.obj.1 ⟶ Z.right.obj.2} {φ' : Z'.right.obj.1 ⟶ Z'.right.obj.2}
    (h : HomPf.mk Z φ = HomPf.mk Z' φ') : repDeg Z φ = repDeg Z' φ' := by
  obtain ⟨V, u, v, huv⟩ := HomPf.eq_iff.mp h
  refine Eq.trans (repDeg_map u φ).symm (Eq.trans ?_ (repDeg_map v φ'))
  exact congrArg (repDeg V) huv

variable {P F} in
theorem repBase_congr {A B : C} {Z Z' : IdxPf P F A B}
    {φ : Z.right.obj.1 ⟶ Z.right.obj.2} {φ' : Z'.right.obj.1 ⟶ Z'.right.obj.2}
    (h : HomPf.mk Z φ = HomPf.mk Z' φ') : repBase Z φ = repBase Z' φ' := by
  obtain ⟨V, u, v, huv⟩ := HomPf.eq_iff.mp h
  refine Eq.trans (repBase_map u φ).symm (Eq.trans ?_ (repBase_map v φ'))
  exact congrArg (repBase V) huv

variable {P F} in
theorem repDiv_congr {A B : C} {Z Z' : IdxPf P F A B}
    {φ : Z.right.obj.1 ⟶ Z.right.obj.2} {φ' : Z'.right.obj.1 ⟶ Z'.right.obj.2}
    (h : HomPf.mk Z φ = HomPf.mk Z' φ') : repDiv Z φ = repDiv Z' φ' := by
  obtain ⟨V, u, v, huv⟩ := HomPf.eq_iff.mp h
  refine Eq.trans (repDiv_map u φ).symm (Eq.trans ?_ (repDiv_map v φ'))
  exact congrArg (repDiv V) huv

/-! ## ★4. `Hom^pf` の上の 3 つの量

★**余錐を使わず、代表元を選んで定義する**(宇宙の制約を避けるため)。 -/

variable {P F} in
/-- ★代表元の選択。 -/
noncomputable def pfRep {A B : C} (f : HomPf P F A B) : IdxPf P F A B :=
  (HomPf.exists_rep f).choose

variable {P F} in
/-- ★選んだ代表元の射。 -/
noncomputable def pfRepMor {A B : C} (f : HomPf P F A B) :
    (pfRep f).right.obj.1 ⟶ (pfRep f).right.obj.2 :=
  (HomPf.exists_rep f).choose_spec.choose

variable {P F} in
theorem pfRep_spec {A B : C} (f : HomPf P F A B) :
    HomPf.mk (pfRep f) (pfRepMor f) = f :=
  (HomPf.exists_rep f).choose_spec.choose_spec

variable {P F} in
/-- ★★**`Hom^pf` の Frobenius 次数**。 -/
noncomputable def pfDeg {A B : C} (f : HomPf P F A B) : ℕ+ := repDeg (pfRep f) (pfRepMor f)

variable {P F} in
/-- ★★**`Hom^pf` の底**。 -/
noncomputable def pfBase {A B : C} (f : HomPf P F A B) :
    (P.toElem.obj A).base ⟶ (P.toElem.obj B).base := repBase (pfRep f) (pfRepMor f)

variable {P F} in
/-- ★★★**`Hom^pf` の零因子**(`Φ^pf` の元)。 -/
noncomputable def pfDiv {A B : C} (f : HomPf P F A B) :
    Pf (Φ.val (P.toElem.obj A).base) := repDiv (pfRep f) (pfRepMor f)

variable {P F} in
@[simp] theorem pfDeg_mk {A B : C} (Z : IdxPf P F A B) (φ : Z.right.obj.1 ⟶ Z.right.obj.2) :
    pfDeg (HomPf.mk Z φ) = repDeg Z φ :=
  repDeg_congr ((pfRep_spec (HomPf.mk Z φ)).trans rfl)

variable {P F} in
@[simp] theorem pfBase_mk {A B : C} (Z : IdxPf P F A B) (φ : Z.right.obj.1 ⟶ Z.right.obj.2) :
    pfBase (HomPf.mk Z φ) = repBase Z φ :=
  repBase_congr ((pfRep_spec (HomPf.mk Z φ)).trans rfl)

variable {P F} in
@[simp] theorem pfDiv_mk {A B : C} (Z : IdxPf P F A B) (φ : Z.right.obj.1 ⟶ Z.right.obj.2) :
    pfDiv (HomPf.mk Z φ) = repDiv Z φ :=
  repDiv_congr ((pfRep_spec (HomPf.mk Z φ)).trans rfl)

/-! ## ★5. 次数と底は合成を保つ

★★**3 つ組の代表元**(`exists_rep3`)で表せば、どちらも `𝔽` の合成規則そのものになる。 -/

variable {P F} in
/-- ★★**次数は合成を保つ**(`𝔽` の規則 `deg(φ≫ψ) = deg ψ · deg φ`)。 -/
theorem pfDeg_comp {A B E : C} (f : HomPf P F A B) (g : HomPf P F B E) :
    pfDeg (compPf P F f g) = pfDeg g * pfDeg f := by
  obtain ⟨V, φ, ψ, rfl, rfl⟩ := exists_rep3 f g
  rw [compPf_mk V φ ψ, pfDeg_mk, pfDeg_mk, pfDeg_mk]
  exact P.degFr_comp φ ψ

/-- ★★**底の合成則**(素の射で述べた形)。 -/
theorem repBaseOf_comp {A B E A' B' E' : C} (a : A ⟶ A') (b : B ⟶ B') (hb : IsIso (P.Base b))
    (e : E ⟶ E') (he : IsIso (P.Base e)) (φ : A' ⟶ B') (ψ : B' ⟶ E') :
    repBaseOf P a e he (φ ≫ ψ) = repBaseOf P a b hb φ ≫ repBaseOf P b e he ψ := by
  haveI := hb
  haveI := he
  show P.Base a ≫ P.Base (φ ≫ ψ) ≫ inv (P.Base e)
    = (P.Base a ≫ P.Base φ ≫ inv (P.Base b)) ≫ (P.Base b ≫ P.Base ψ ≫ inv (P.Base e))
  simp only [P.Base_comp, Category.assoc, IsIso.inv_hom_id_assoc]

variable {P F} in
/-- ★★**底は合成を保つ**。 -/
theorem pfBase_comp {A B E : C} (f : HomPf P F A B) (g : HomPf P F B E) :
    pfBase (compPf P F f g) = pfBase f ≫ pfBase g := by
  obtain ⟨V, φ, ψ, rfl, rfl⟩ := exists_rep3 f g
  rw [compPf_mk V φ ψ, pfBase_mk, pfBase_mk, pfBase_mk]
  exact repBaseOf_comp P V.hom.hom.1 V.hom.hom.2.1
    ((idx12 P F A B E).obj V).hom.property.2.1.2 V.hom.hom.2.2
    ((idx23 P F A B E).obj V).hom.property.2.1.2 φ ψ

/-! ## ★6. 零因子の合成則

★★`𝔽` の規則 `Div(φ≫ψ) = (Base φ)^*(Div ψ) + deg(ψ)·Div φ` の **`Pf` 版**。
★**分母(根の次数)が 3 脚で共通**であることが効く。

★★**実装上の抜け道 2 つ(2026-08-17 の測定)**:
1. 目標の**内側**を `rw` / `simp only` で書き換えようとすると `Φ.map` の
   coercion の形が食い違って当たらない。★**`Pf.mk` の分子と分母を
   `congrArg₂ Pf.mk` でまとめて `have` にし、外側で `rw` する**と通る。
2. ★**分母共通版 `Pf.mk_add_mk_same`(★0)を先に置く** —— これで
   `Pf.sound` を目標の中で使わずに済む。 -/

variable {P F} in
/-- ★★★**零因子は合成を保つ**(`𝔽_{Φ^pf}` の合成規則)。 -/
theorem pfDiv_comp {A B E : C} (f : HomPf P F A B) (g : HomPf P F B E) :
    pfDiv (compPf P F f g)
      = Pf.map (Φ.map (pfBase f)) (pfDiv g) + ((pfDeg g : ℕ+) : ℕ) • pfDiv f := by
  obtain ⟨V, φ, ψ, rfl, rfl⟩ := exists_rep3 f g
  haveI hb : IsIso (P.Base V.hom.hom.2.1) := ((idx12 P F A B E).obj V).hom.property.2.1.2
  rw [compPf_mk V φ ψ, pfDiv_mk, pfDiv_mk, pfDiv_mk, pfBase_mk, pfDeg_mk]
  have hmn : P.degFr ((idx23 P F A B E).obj V).hom.hom.1
      = P.degFr ((idx13 P F A B E).obj V).hom.hom.1 := V.hom.property.2.2.2.1.symm
  have hcomp : repBase ((idx12 P F A B E).obj V) φ
        ≫ P.Base ((idx23 P F A B E).obj V).hom.hom.1
      = P.Base ((idx12 P F A B E).obj V).hom.hom.1 ≫ P.Base φ := by
    show (P.Base V.hom.hom.1 ≫ P.Base φ ≫ inv (P.Base V.hom.hom.2.1))
        ≫ P.Base V.hom.hom.2.1 = _
    simp only [Category.assoc, IsIso.inv_hom_id, Category.comp_id]
    rfl
  have hY : Φ.map (repBase ((idx12 P F A B E).obj V) φ)
      (Φ.map (P.Base ((idx23 P F A B E).obj V).hom.hom.1) (P.Div ψ))
      = Φ.map (P.Base ((idx13 P F A B E).obj V).hom.hom.1)
        (Φ.map (P.Base φ) (P.Div ψ)) :=
    ((Φ.map_comp (P.Base ((idx23 P F A B E).obj V).hom.hom.1)
        (repBase ((idx12 P F A B E).obj V) φ) (P.Div ψ)).symm.trans
      (congrArg (fun t => Φ.map t (P.Div ψ)) hcomp)).trans
      (Φ.map_comp (P.Base φ) (P.Base ((idx13 P F A B E).obj V).hom.hom.1) (P.Div ψ))
  have hPf : Pf.map (Φ.map (repBase ((idx12 P F A B E).obj V) φ))
        (repDiv ((idx23 P F A B E).obj V) ψ)
      = Pf.mk (Φ.map (P.Base ((idx13 P F A B E).obj V).hom.hom.1)
          (Φ.map (P.Base φ) (P.Div ψ))) (repRoot ((idx13 P F A B E).obj V)) := by
    show Pf.map _ (Pf.mk _ _) = _
    rw [Pf.map_mk]
    exact congrArg₂ Pf.mk hY hmn
  have hSm : ((repDeg ((idx23 P F A B E).obj V) ψ : ℕ+) : ℕ)
        • repDiv ((idx12 P F A B E).obj V) φ
      = Pf.mk (((P.degFr ψ : ℕ+) : ℕ)
          • Φ.map (P.Base ((idx13 P F A B E).obj V).hom.hom.1) (P.Div φ))
        (repRoot ((idx13 P F A B E).obj V)) := by
    show _ • Pf.mk _ _ = _
    rw [Pf.nsmul_mk]
    rfl
  rw [hPf, hSm]
  refine Eq.trans ?_ (Pf.mk_add_mk_same
    (Φ.map (P.Base ((idx13 P F A B E).obj V).hom.hom.1) (Φ.map (P.Base φ) (P.Div ψ)))
    (((P.degFr ψ : ℕ+) : ℕ)
      • Φ.map (P.Base ((idx13 P F A B E).obj V).hom.hom.1) (P.Div φ))
    (repRoot ((idx13 P F A B E).obj V))).symm
  simp only [repDiv, repRoot]
  rw [P.Div_comp, map_add, map_nsmul]
  rfl

/-! ## ★6-b. `𝒞 → 𝒞^pf` の像での計算則

★`toHomPf φ` は添字 `idxOne`(両脚とも `𝟙`)での代表元なので、
3 つの量は**そのまま `𝒞` の値**になる。★零因子だけは分母 `1` が付く。 -/

variable {P F} in
@[simp] theorem pfDeg_toHomPf {A B : C} (φ : A ⟶ B) :
    pfDeg (toHomPf (F := F) φ) = P.degFr φ := by
  rw [toHomPf, pfDeg_mk]
  rfl

variable {P F} in
@[simp] theorem pfBase_toHomPf {A B : C} (φ : A ⟶ B) :
    pfBase (toHomPf (F := F) φ) = P.Base φ := by
  haveI : IsIso (P.Base (𝟙 B)) := by rw [P.Base_id]; infer_instance
  rw [toHomPf, pfBase_mk]
  show P.Base (𝟙 A) ≫ P.Base φ ≫ inv (P.Base (𝟙 B)) = P.Base φ
  rw [← Category.assoc, IsIso.comp_inv_eq, P.Base_id, P.Base_id]
  simp

variable {P F} in
@[simp] theorem pfDiv_toHomPf {A B : C} (φ : A ⟶ B) :
    pfDiv (toHomPf (F := F) φ) = Pf.mk (P.Div φ) 1 := by
  rw [toHomPf, pfDiv_mk]
  show Pf.mk (Φ.map (P.Base (𝟙 A)) (P.Div φ)) (P.degFr (𝟙 A)) = _
  rw [P.Base_id, P.degFr_id]
  exact congrArg (fun t => Pf.mk t 1) (Φ.map_id _ _)

/-! ## ★6-c. ★★**関手 `𝒞^pf ⥤ 𝔽_{Φ^pf}`**

原文 (FrdI p.59):
> the left; the lower horizontal arrow is induced by the natural morphism of monoids

★★**`Φ` が divisorial なので `Φ.val A` は sharp**、したがって `Φ^pf` が作れる。
★対象は `𝒟` 側で変わらない(`(A,n)` の底は `A` の底)。
★射は `⟨pfBase f, pfDiv f, pfDeg f⟩`——★**関手則は ★5・★6 の 3 本がそのまま与える**。

★★ここでは **`PfCat`(`n = 1` の部分)版**を組む。
原文の `𝒞^pf` は対象が対 `(A,n)` の `pfRootCategory` なので、
★**これはその「`n=1` への制限」である**(★7-b' を見よ)。 -/

include P in
/-- ★`Φ` が divisorial なら、各 `Φ.val A` は sharp。 -/
theorem phiSharp (A : D) : IsSharp (Φ.val A) := (P.divisorial A).2

/-- ★★★**[FrdI] Proposition 3.2, (i)** —— `𝒞^pf → 𝔽_{Φ^pf}`(`n = 1` 版)。 -/
noncomputable def pfToElem : PfCat P F ⥤ ElemFrobCat (Φ.pfOn (phiSharp P)) where
  obj A := ⟨(P.toElem.obj (pfDown P F A)).base⟩
  map f := ⟨pfBase f, pfDiv f, pfDeg f⟩
  map_id A := by
    refine ElemFrobCat.Hom.ext ?_ ?_ ?_
    · show pfBase (toHomPf (F := F) (𝟙 (pfDown P F A))) = _
      rw [pfBase_toHomPf, P.Base_id]
      rfl
    · show pfDiv (toHomPf (F := F) (𝟙 (pfDown P F A))) = _
      rw [pfDiv_toHomPf, P.Div_id]
      rfl
    · show pfDeg (toHomPf (F := F) (𝟙 (pfDown P F A))) = _
      rw [pfDeg_toHomPf, P.degFr_id]
      rfl
  map_comp f g := by
    refine ElemFrobCat.Hom.ext ?_ ?_ ?_
    · exact pfBase_comp f g
    · exact pfDiv_comp f g
    · exact pfDeg_comp f g

/-! ## ★6-d. ★★★**1-可換図式**([FrdI] Proposition 3.2, (i))

原文 (FrdI p.59):
> the left; the lower horizontal arrow is induced by the natural morphism of monoids

★4 つの関手:

| 辺 | 関手 |
|---|---|
| 上(`𝒞 → 𝒞^pf`) | `toPfCat`(`Definition 3.1, (iii)`) |
| 左(`𝒞 → 𝔽_Φ`) | `P.toElem`(pre-Frobenioid 構造そのもの) |
| 右(`𝒞^pf → 𝔽_{Φ^pf}`) | `pfToElem`(★6-c) |
| 下(`𝔽_Φ → 𝔽_{Φ^pf}`) | `elemFrobMap` を `Φ → Φ^pf` に当てたもの |

★★**下の辺の「自然なモノイドの射 `Φ → Φ^pf`」は `Pf.of`(`m ↦ m/1`)**である。 -/

/-- ★★**自然な射 `Φ ⟶ Φ^pf`** —— 各点で `Pf.of`(`m ↦ m/1`)。

★自然性は `Pf.map f (m/1) = (f m)/1`(`Pf.map_mk`)そのものである。 -/
noncomputable def phiToPf : Φ.functor ⟶ (Φ.pfOn (phiSharp P)).functor where
  app X := AddCommMonCat.ofHom (Pf.of (M := Φ.val X.unop))
  naturality X Y f := by
    refine AddCommMonCat.ext (fun m => ?_)
    show Pf.of (Φ.map f.unop m) = Pf.map (Φ.map f.unop) (Pf.of m)
    rw [Pf.of_apply, Pf.of_apply, Pf.map_mk]

/-- ★★**下の辺** `𝔽_Φ ⥤ 𝔽_{Φ^pf}`。 -/
noncomputable def elemToPfElem : ElemFrobCat Φ ⥤ ElemFrobCat (Φ.pfOn (phiSharp P)) :=
  ElemFrobCat.elemFrobMap (phiToPf P)

/-- ★★★**[FrdI] Proposition 3.2, (i)** —— **図式は 1-可換**である。

★★実際には**厳密に等しい**(`Iso.refl` で足りる)ことが分かった:
対象は両側とも `⟨(P.toElem.obj A).base⟩`、射は両側とも
`⟨P.Base φ, Pf.mk (P.Div φ) 1, P.degFr φ⟩` になる。 -/
theorem pfSquare_obj (A : C) :
    (toPfCat P F ⋙ pfToElem P F).obj A = (P.toElem ⋙ elemToPfElem P).obj A := rfl

/-- ★★★**[FrdI] Proposition 3.2, (i)** —— 射の側の 1-可換性。 -/
theorem pfSquare_map {A B : C} (φ : A ⟶ B) :
    (toPfCat P F ⋙ pfToElem P F).map φ = (P.toElem ⋙ elemToPfElem P).map φ := by
  refine ElemFrobCat.Hom.ext ?_ ?_ ?_
  · exact pfBase_toHomPf (F := F) φ
  · exact pfDiv_toHomPf (F := F) φ
  · exact pfDeg_toHomPf (F := F) φ

/-- ★★★**[FrdI] Proposition 3.2, (i)** —— 図式の 1-可換性(関手の同型として)。 -/
noncomputable def pfSquare : toPfCat P F ⋙ pfToElem P F ≅ P.toElem ⋙ elemToPfElem P :=
  NatIso.ofComponents (fun X => Iso.refl ((toPfCat P F ⋙ pfToElem P F).obj X)) (by
    intro A B φ
    show (toPfCat P F ⋙ pfToElem P F).map φ ≫ 𝟙 _
      = 𝟙 _ ≫ (P.toElem ⋙ elemToPfElem P).map φ
    rw [Category.comp_id, Category.id_comp]
    exact pfSquare_map P F φ)

/-! ## ★6-e. 段 D の残り —— `𝒞^pf` の連結性

★★**`PfCat P F` は定義からして `C` そのもの**(`def PfCat _P _F : Type u2 := C`)で、
`toPfCat` は**対象について恒等**である。★したがって `C` の zigzag を
`toPfCat` で送るだけでよい。 -/

/-- ★**`𝒞^pf` は連結**(`n = 1` 版)。 -/
theorem pfCat_isConnected : IsConnected (PfCat P F) := by
  haveI := P.connectedC
  haveI : Nonempty (PfCat P F) := (inferInstance : Nonempty C)
  refine zigzag_isConnected ?_
  intro A B
  exact zigzag_obj_of_zigzag (toPfCat P F)
    (isPreconnected_zigzag (J := C) (pfDown P F A) (pfDown P F B))

/-! ## ★6-f. 段 D の最後 —— `𝒞^pf` が totally epimorphic

★`IsTotallyEpimorphic` は「**すべての射が epi**」。
★★**筋**: `f` `g` `h` を 1 つの 3 脚添字の上に揃え、`compPf_mk` で
`mk (idx13 V) (φ≫ψ) = mk (idx13 V) (φ≫χ)` に落とし、
`idxTransport_comp_pair` で遷移を分けてから ★**`P.totEpiC` で `φ` の遷移像を消す**。

★★**添字圏が細い(`idx_hom_ext` / `idx3_hom_ext`)ことが 2 度効く** ——
`HomPf.eq_iff` が返す 2 本の遷移射を同一視するときと、
`idx13` の像に持ち上げた遷移射を `idx13.map t` と同一視するときである。 -/

/-- ★**cofinal な関手からは、任意の対象へ「上に行く射」が取れる**。 -/
theorem exists_hom_of_final {J K : Type*} [Category J] [Category K]
    (G : J ⥤ K) [G.Final] (Z : K) : ∃ V : J, Nonempty (Z ⟶ G.obj V) := by
  obtain ⟨s⟩ := (Functor.Final.out (F := G) Z).is_nonempty
  exact ⟨s.right, ⟨s.hom⟩⟩

variable {P F} in
/-- ★★**平行 2 射も同じ 3 脚添字の上に揃う** —— `exists_rep3` を 2 回使い、
`IdxPf3` の filtered 性で上界を取るだけ。 -/
theorem exists_rep3_pair {A B E : C} (f : HomPf P F A B) (g h : HomPf P F B E) :
    ∃ (V : IdxPf3 P F A B E) (φ : V.right.obj.1 ⟶ V.right.obj.2.1)
      (ψ χ : V.right.obj.2.1 ⟶ V.right.obj.2.2),
      f = HomPf.mk ((idx12 P F A B E).obj V) φ ∧
      g = HomPf.mk ((idx23 P F A B E).obj V) ψ ∧
      h = HomPf.mk ((idx23 P F A B E).obj V) χ := by
  obtain ⟨V₁, φ₁, ψ₁, hf₁, hg⟩ := exists_rep3 f g
  obtain ⟨V₂, _φ₂, χ₂, _hf₂, hh⟩ := exists_rep3 f h
  refine ⟨IsFiltered.max V₁ V₂,
    idxTransport P F ((idx12 P F A B E).map (IsFiltered.leftToMax V₁ V₂)) φ₁,
    idxTransport P F ((idx23 P F A B E).map (IsFiltered.leftToMax V₁ V₂)) ψ₁,
    idxTransport P F ((idx23 P F A B E).map (IsFiltered.rightToMax V₁ V₂)) χ₂,
    ?_, ?_, ?_⟩
  · rw [HomPf.mk_map]; exact hf₁
  · rw [HomPf.mk_map]; exact hg
  · rw [HomPf.mk_map]; exact hh

variable {P F} in
/-- ★★**`Hom^pf` の左簡約** —— これが `𝒞^pf` の全射性の本体。 -/
theorem homPf_epi_cancel {A B E : C} (f : HomPf P F A B) (g h : HomPf P F B E)
    (hgh : compPf P F f g = compPf P F f h) : g = h := by
  obtain ⟨V, φ, ψ, χ, rfl, rfl, rfl⟩ := exists_rep3_pair f g h
  rw [compPf_mk V φ ψ, compPf_mk V φ χ] at hgh
  obtain ⟨W, u, v, hw⟩ := HomPf.eq_iff.mp hgh
  rw [idx_hom_ext u v] at hw
  -- `W` を `idx13` の像へ持ち上げ、さらに `V` と共通の上界を取る
  obtain ⟨V₀, ⟨k⟩⟩ := exists_hom_of_final (idx13 P F A B E) W
  set V' := IsFiltered.max V V₀ with hV'
  set t : V ⟶ V' := IsFiltered.leftToMax V V₀ with ht
  set t₀ : V₀ ⟶ V' := IsFiltered.rightToMax V V₀ with ht₀
  have hsplit : (idx13 P F A B E).map t
      = v ≫ k ≫ (idx13 P F A B E).map t₀ := idx_hom_ext _ _
  have key : idxTransport P F ((idx13 P F A B E).map t) (φ ≫ ψ)
      = idxTransport P F ((idx13 P F A B E).map t) (φ ≫ χ) := by
    rw [hsplit]
    simp only [idxTransport_comp]
    rw [hw]
  rw [← idxTransport_comp_pair t φ ψ, ← idxTransport_comp_pair t φ χ] at key
  haveI : Epi (idxTransport P F ((idx12 P F A B E).map t) φ) := P.totEpiC _ _ _
  have hfin := (cancel_epi (idxTransport P F ((idx12 P F A B E).map t) φ)).mp key
  rw [← HomPf.mk_map ((idx23 P F A B E).map t) ψ, ← HomPf.mk_map ((idx23 P F A B E).map t) χ,
    hfin]

/-- ★★**`𝒞^pf` は totally epimorphic**(`n = 1` 版)。 -/
theorem pfCat_totEpi : IsTotallyEpimorphic (PfCat P F) :=
  fun _ _ f => ⟨fun g h hgh => homPf_epi_cancel f g h hgh⟩

/-! ## ★6-g. ★★★**`𝒞^pf` の pre-Frobenioid 構造**

原文 (FrdI p.59):
> the left; the lower horizontal arrow is induced by the natural morphism of monoids

★★原文の「In particular, the functor Cpf to FPhipf determines a pre-Frobenioid
structure on Cpf」の実体。★**6 フィールドの内訳**:

| # | フィールド | 出どころ |
|---|---|---|
| 1 | `toElem` | `pfToElem`(★6-c) |
| 2 | `divisorial` | `Pf.isDivisorial'`(★0-b) |
| 3 | `totEpiC` | `pfCat_totEpi`(★6-f) |
| 4 | `totEpiD` | ★**`P.totEpiD` そのまま**(`𝒟` は完全化で変わらない) |
| 5 | `connectedC` | `pfCat_isConnected`(★6-e) |
| 6 | `connectedD` | ★**`P.connectedD` そのまま** |
-/

/-- ★★★**[FrdI] Proposition 3.2, (i)** —— `𝒞^pf` の pre-Frobenioid 構造(`n = 1` 版)。 -/
noncomputable def pfPre : PreFrobenioid (PfCat P F) (Φ.pfOn (phiSharp P)) where
  toElem := pfToElem P F
  divisorial A := Pf.isDivisorial' (P.divisorial A)
  totEpiC := pfCat_totEpi P F
  totEpiD := P.totEpiD
  connectedC := pfCat_isConnected P F
  connectedD := P.connectedD

/-! ## ★6-h. 段 B′ の下ごしらえ —— **根の取り替え(`rootIso`)で 3 量はどう動くか**

★★原文の `𝒞^pf` は対象が対 `(A,n)` の圏(`pfRootCategory`)であり、
その射は `Hom^pf(A^{(m)}, B^{(n)})` である。★合成 `compRoot` は
`rtRootIso`(= `rootIso`)で根を揃えてから `compPf` するので、
★**`rootIso` の下での `pfDeg` / `pfBase` / `pfDiv` の変化則**が要る。

★`rootIso_hom_mk` は「**代表元の射そのものは変わらず、添字だけ押し出される**」と言う。
押し出された添字の脚は `a ≫ (元の脚)`、`b ≫ (元の脚)` である。したがって:

| 量 | 変化 |
|---|---|
| `pfDeg` | ★★**不変**(`repDeg Z φ = P.degFr φ` は添字に依らない) |
| `pfBase` | `Base a ≫ (元) ≫ inv (Base b)` で共役 |
| `pfDiv` | 分子は `Φ.map (Base a)` で押し出し、★**分母は `degFr a` 倍**になる |

★★**分母が増えるのが要点**である —— `(A,n)` は「`A` の `n` 乗根」なので、
零因子も `n` で割られる。 -/

variable {P F} in
/-- ★★**根を取り替えても Frobenius 次数は変わらない**。 -/
theorem pfDeg_rootIso {A B A' B' : C} (a : A ⟶ A') (ha : IsFrobeniusType P a)
    (b : B ⟶ B') (hb : IsFrobeniusType P b) (hd : P.degFr a = P.degFr b)
    (z : HomPf P F A' B') :
    pfDeg ((rootIso (F := F) a ha b hb hd).hom z) = pfDeg z := by
  obtain ⟨V, φ, rfl⟩ := HomPf.exists_rep z
  rw [rootIso_hom_mk, pfDeg_mk, pfDeg_mk]
  rfl

variable {P F} in
/-- ★**押し出した添字の根の次数**は `degFr a` 倍になる。 -/
theorem repRoot_pushIdx {A B A' B' : C} (a : A ⟶ A') (ha : IsFrobeniusType P a)
    (b : B ⟶ B') (hb : IsFrobeniusType P b) (hd : P.degFr a = P.degFr b)
    (V : IdxPf P F A' B') :
    repRoot ((pushIdx (F := F) a ha b hb hd).obj V) = repRoot V * P.degFr a :=
  P.degFr_comp a V.hom.hom.1

variable {P F} in
/-- ★**押し出した添字の底**は `Base a` と `(Base b)⁻¹` で共役を取ったもの。 -/
theorem repBase_pushIdx {A B A' B' : C} (a : A ⟶ A') (ha : IsFrobeniusType P a)
    (b : B ⟶ B') (hb : IsFrobeniusType P b) (hd : P.degFr a = P.degFr b)
    (V : IdxPf P F A' B') (φ : V.right.obj.1 ⟶ V.right.obj.2) :
    repBase ((pushIdx (F := F) a ha b hb hd).obj V) φ ≫ P.Base b
      = P.Base a ≫ repBase V φ := by
  haveI hp : IsIso (P.Base (b ≫ V.hom.hom.2)) :=
    ((pushIdx (F := F) a ha b hb hd).obj V).hom.property.2.1.2
  haveI h2 : IsIso (P.Base V.hom.hom.2) := V.hom.property.2.1.2
  refine (cancel_mono (P.Base V.hom.hom.2)).mp ?_
  show ((P.Base (a ≫ V.hom.hom.1) ≫ P.Base φ ≫ inv (P.Base (b ≫ V.hom.hom.2)))
        ≫ P.Base b) ≫ P.Base V.hom.hom.2
      = (P.Base a ≫ (P.Base V.hom.hom.1 ≫ P.Base φ ≫ inv (P.Base V.hom.hom.2)))
        ≫ P.Base V.hom.hom.2
  simp only [Category.assoc]
  rw [← P.Base_comp b V.hom.hom.2, IsIso.inv_hom_id, Category.comp_id,
    IsIso.inv_hom_id, Category.comp_id, P.Base_comp]
  simp only [Category.assoc]

variable {P F} in
/-- ★★**押し出した添字の零因子** —— 分子は `Φ.map (Base a)` で押し出され、
★**分母は `degFr a` 倍**になる。★これが「`(A,n)` は `A` の `n` 乗根」の正体。 -/
theorem repDiv_pushIdx {A B A' B' : C} (a : A ⟶ A') (ha : IsFrobeniusType P a)
    (b : B ⟶ B') (hb : IsFrobeniusType P b) (hd : P.degFr a = P.degFr b)
    (V : IdxPf P F A' B') (φ : V.right.obj.1 ⟶ V.right.obj.2) :
    repDiv ((pushIdx (F := F) a ha b hb hd).obj V) φ
      = Pf.mk (Φ.map (P.Base a) (Φ.map (P.Base V.hom.hom.1) (P.Div φ)))
        (repRoot V * P.degFr a) := by
  show Pf.mk (Φ.map (P.Base (a ≫ V.hom.hom.1)) (P.Div φ)) (P.degFr (a ≫ V.hom.hom.1)) = _
  rw [P.Base_comp]
  exact congrArg₂ Pf.mk (Φ.map_comp (P.Base V.hom.hom.1) (P.Base a) (P.Div φ))
    (P.degFr_comp a V.hom.hom.1)

variable {P F} in
/-- ★★**根を取り替えると底は共役される**。 -/
theorem pfBase_rootIso {A B A' B' : C} (a : A ⟶ A') (ha : IsFrobeniusType P a)
    (b : B ⟶ B') (hb : IsFrobeniusType P b) (hd : P.degFr a = P.degFr b)
    (z : HomPf P F A' B') :
    pfBase ((rootIso (F := F) a ha b hb hd).hom z) ≫ P.Base b
      = P.Base a ≫ pfBase z := by
  obtain ⟨V, φ, rfl⟩ := HomPf.exists_rep z
  rw [rootIso_hom_mk, pfBase_mk, pfBase_mk]
  exact repBase_pushIdx a ha b hb hd V φ

variable {P F} in
/-- ★★★**根を取り替えると零因子は押し出されたうえで `degFr a` で割られる**。

★★これが `Definition 3.1, (iii)` の「`(A,n)` は `A` の `n` 乗根と思え」の
零因子側の内容である。 -/
theorem pfDiv_rootIso {A B A' B' : C} (a : A ⟶ A') (ha : IsFrobeniusType P a)
    (b : B ⟶ B') (hb : IsFrobeniusType P b) (hd : P.degFr a = P.degFr b)
    (z : HomPf P F A' B') :
    pfDiv ((rootIso (F := F) a ha b hb hd).hom z)
      = Pf.divBy (P.degFr a) (Pf.map (Φ.map (P.Base a)) (pfDiv z)) := by
  obtain ⟨V, φ, rfl⟩ := HomPf.exists_rep z
  rw [rootIso_hom_mk, pfDiv_mk, pfDiv_mk, repDiv_pushIdx]
  show _ = Pf.divBy (P.degFr a) (Pf.map (Φ.map (P.Base a)) (Pf.mk _ _))
  rw [Pf.map_mk, Pf.divBy_mk]
  exact congrArg (Pf.mk _) (mul_comm _ _)

/-! ## ★6-i. 段 B′ —— 対象が対 `(A,n)` の `𝒞^pf`(`pfRootCategory`)

原文 (FrdI p.57):
> follows: The objects of Cpf are pairs (A, n), where A

★★`HomRoot (A,n) (B,m) = Hom^pf(A^{(m)}, B^{(n)})` である。
★3 量の定義は「**`(A,n)` 側へ引き戻す**」形になる:

| 量 | 定義 | 不変性の根拠 |
|---|---|---|
| `rootDeg` | `pfDeg f` | ★`pfDeg` は根の取り替えで不変 |
| `rootBase` | `Base(A→A^{(m)}) ≫ pfBase f ≫ Base(B→B^{(n)})⁻¹` | 底同型の合成 |
| `rootDiv` | ★**`(押し出した pfDiv) / (n·m)`** | ★★合成則が正規化を一意に決める |

★★**`rootDiv` で割る数は `n·m = X.root * Y.root`** である。
★★**これは合成則から逆算して確定した**(下の `rootDiv` の docstring を見よ)。 -/

variable {P F} in
/-- ★根の取り替え(選んだ版)で次数は不変。 -/
theorem pfDeg_rtRootIso_hom (A B : C) {dA dB e tA tB : ℕ+} (hA : tA = e * dA) (hB : tB = e * dB)
    (z : HomPf P F (rtObj P F A tA) (rtObj P F B tB)) :
    pfDeg ((rtRootIso P F A B hA hB).hom z) = pfDeg z :=
  pfDeg_rootIso _ _ _ _ _ z

variable {P F} in
/-- ★逆向きも同じ。 -/
theorem pfDeg_rtRootIso_inv (A B : C) {dA dB e tA tB : ℕ+} (hA : tA = e * dA) (hB : tB = e * dB)
    (z : HomPf P F (rtObj P F A dA) (rtObj P F B dB)) :
    pfDeg ((rtRootIso P F A B hA hB).inv z) = pfDeg z := by
  have h := pfDeg_rtRootIso_hom A B hA hB ((rtRootIso P F A B hA hB).inv z)
  rw [Iso.inv_hom_id_apply] at h
  exact h.symm

variable {P F} in
/-- ★★**`(A,n) ⟶ (B,m)` の Frobenius 次数**。 -/
noncomputable def rootDeg {X Y : PfRootObj P F} (f : HomRoot P F X Y) : ℕ+ := pfDeg f

variable {P F} in
/-- ★★**`(A,n) ⟶ (B,m)` の底** —— `A^{(m)}` 側の底を `A` 側へ共役で戻す。 -/
noncomputable def rootBase {X Y : PfRootObj P F} (f : HomRoot P F X Y) :
    (P.toElem.obj X.obj).base ⟶ (P.toElem.obj Y.obj).base :=
  letI : IsIso (P.Base (rtExt P F Y.obj X.root)) := (rtExt_frobType P F Y.obj X.root).2
  P.Base (rtExt P F X.obj Y.root) ≫ pfBase f ≫ inv (P.Base (rtExt P F Y.obj X.root))

variable {P F} in
/-- ★★★**`(A,n) ⟶ (B,m)` の零因子** —— 押し出してから ★**`n·m` で割る**。

★★**割る数が `n·m` である理由**(2026-08-17、合成則から逆算して確定):
`(A,n)` は `A^{1/n}`、`(B,m)` は `B^{1/m}` と思うので、
代表の `A^{(m)} ⟶ B^{(n)}` は求める射の **`n·m` 乗**にあたる。
★したがって零因子は `n·m` で割られる。

★★**一度 `m` だけで割る形で書いたが、それでは合成則が成り立たない**
(第 1 項の係数が左辺 `n`・右辺 `m` になって食い違う)。★**合成則が正しい正規化を
一意に決める**——`c_{X,Y} · Y.root = c_{Y,Z} · X.root` と
`c_{X,Y} · Y.root = c_{X,Y} · Z.root` の 2 条件から `c_{X,Y} = X.root · Y.root`。 -/
noncomputable def rootDiv {X Y : PfRootObj P F} (f : HomRoot P F X Y) :
    Pf (Φ.val (P.toElem.obj X.obj).base) :=
  Pf.divBy (X.root * Y.root)
    (Pf.map (Φ.map (P.Base (rtExt P F X.obj Y.root))) (pfDiv f))

variable {P F} in
/-- ★★**根つきでも次数は合成を保つ** —— `compRoot` の 3 つの `rtRootIso` は
どれも次数を動かさないので、`pfDeg_comp` がそのまま効く。 -/
theorem rootDeg_comp {X Y Z : PfRootObj P F} (f : HomRoot P F X Y) (g : HomRoot P F Y Z) :
    rootDeg (compRoot P F f g) = rootDeg g * rootDeg f := by
  show pfDeg (compRoot P F f g) = _
  rw [compRoot, pfDeg_rtRootIso_hom, pfDeg_comp, pfDeg_rtRootIso_inv, pfDeg_rtRootIso_inv]
  rfl

variable {P F} in
/-- ★根の取り替え(選んだ版)での底の変化 —— **`inv` を使わない形**。 -/
theorem pfBase_rtRootIso (A B : C) {dA dB e tA tB : ℕ+} (hA : tA = e * dA) (hB : tB = e * dB)
    (w : HomPf P F (rtObj P F A tA) (rtObj P F B tB)) :
    pfBase ((rtRootIso P F A B hA hB).hom w) ≫ P.Base (rtLift P F B hB)
      = P.Base (rtLift P F A hA) ≫ pfBase w :=
  pfBase_rootIso _ _ _ _ _ w

variable {P F} in
/-- ★`rootBase` の**特徴づけ**(`inv` を使わない形)。 -/
theorem rootBase_spec {X Y : PfRootObj P F} (f : HomRoot P F X Y) :
    rootBase f ≫ P.Base (rtExt P F Y.obj X.root)
      = P.Base (rtExt P F X.obj Y.root) ≫ pfBase f := by
  haveI : IsIso (P.Base (rtExt P F Y.obj X.root)) := (rtExt_frobType P F Y.obj X.root).2
  show (P.Base (rtExt P F X.obj Y.root) ≫ pfBase f
      ≫ inv (P.Base (rtExt P F Y.obj X.root))) ≫ P.Base (rtExt P F Y.obj X.root) = _
  rw [Category.assoc, Category.assoc, IsIso.inv_hom_id, Category.comp_id]

variable {P F} in
/-- ★`rootBase` は上の特徴づけで**一意に決まる**(`Base(eB)` が同型ゆえ mono)。 -/
theorem rootBase_uniq {X Y : PfRootObj P F} (f : HomRoot P F X Y)
    (u : (P.toElem.obj X.obj).base ⟶ (P.toElem.obj Y.obj).base)
    (h : u ≫ P.Base (rtExt P F Y.obj X.root)
      = P.Base (rtExt P F X.obj Y.root) ≫ pfBase f) : u = rootBase f := by
  haveI : IsIso (P.Base (rtExt P F Y.obj X.root)) := (rtExt_frobType P F Y.obj X.root).2
  exact (cancel_mono (P.Base (rtExt P F Y.obj X.root))).mp (h.trans (rootBase_spec f).symm)

variable {P F} in
/-- ★根を**上げる**向きの底の変化。 -/
theorem pfBase_rtRootIso_inv (A B : C) {dA dB e tA tB : ℕ+} (hA : tA = e * dA)
    (hB : tB = e * dB) (z : HomPf P F (rtObj P F A dA) (rtObj P F B dB)) :
    pfBase z ≫ P.Base (rtLift P F B hB)
      = P.Base (rtLift P F A hA) ≫ pfBase ((rtRootIso P F A B hA hB).inv z) := by
  have h := pfBase_rtRootIso A B hA hB ((rtRootIso P F A B hA hB).inv z)
  rw [Iso.inv_hom_id_apply] at h
  exact h

variable {P F} in
/-- ★★**根つきでも底は合成を保つ**。 -/
theorem rootBase_comp {X Y Z : PfRootObj P F} (f : HomRoot P F X Y) (g : HomRoot P F Y Z) :
    rootBase (compRoot P F f g) = rootBase f ≫ rootBase g := by
  refine (rootBase_uniq (compRoot P F f g) (rootBase f ≫ rootBase g) ?_).symm
  haveI : IsIso (P.Base (rtLift P F Z.obj
      (show Y.root * X.root = Y.root * X.root from rfl))) :=
    (rtLift_frobType P F Z.obj _).2
  refine (cancel_mono (P.Base (rtLift P F Z.obj
      (show Y.root * X.root = Y.root * X.root from rfl)))).mp ?_
  simp only [compRoot]
  conv_rhs =>
    rw [Category.assoc, pfBase_rtRootIso, ← Category.assoc, ← P.Base_comp, rtLift_ext,
      pfBase_comp]
  have hg : rootBase g ≫ P.Base (rtExt P F Z.obj (Y.root * X.root))
      = P.Base (rtExt P F Y.obj (Z.root * X.root))
        ≫ pfBase ((rtRootIso P F Y.obj Z.obj
            (show Z.root * X.root = X.root * Z.root from mul_comm _ _)
            (show Y.root * X.root = X.root * Y.root from mul_comm _ _)).inv g) := by
    rw [← rtLift_ext P F Z.obj (show Y.root * X.root = X.root * Y.root from mul_comm _ _),
      P.Base_comp, ← Category.assoc, rootBase_spec, Category.assoc,
      pfBase_rtRootIso_inv, ← Category.assoc, ← P.Base_comp, rtLift_ext]
  have hf : rootBase f ≫ P.Base (rtExt P F Y.obj (Z.root * X.root))
      = P.Base (rtExt P F X.obj (Z.root * Y.root))
        ≫ pfBase ((rtRootIso P F X.obj Y.obj
            (show Z.root * Y.root = Z.root * Y.root from rfl)
            (show Z.root * X.root = Z.root * X.root from rfl)).inv f) := by
    rw [← rtLift_ext P F Y.obj (show Z.root * X.root = Z.root * X.root from rfl),
      P.Base_comp, ← Category.assoc, rootBase_spec, Category.assoc,
      pfBase_rtRootIso_inv, ← Category.assoc, ← P.Base_comp, rtLift_ext]
  conv_lhs =>
    rw [Category.assoc, Category.assoc, ← P.Base_comp, rtLift_ext, hg, ← Category.assoc,
      hf, Category.assoc]

/-- ★`Pf.map` は合成を保つ。 -/
theorem Pf.map_map {M N K : Type w} [AddCommMonoid M] [AddCommMonoid N] [AddCommMonoid K]
    (f : N →+ K) (g : M →+ N) (x : Pf M) : Pf.map f (Pf.map g x) = Pf.map (f.comp g) x := by
  induction x using Pf.inductionOn with | _ m a => rfl

variable {P F} in
/-- ★`rootDiv` の**特徴づけ**(`divBy` を掛け算に直した形)。 -/
theorem rootDiv_spec {X Y : PfRootObj P F} (f : HomRoot P F X Y) :
    ((X.root * Y.root : ℕ+) : ℕ) • rootDiv f
      = Pf.map (Φ.map (P.Base (rtExt P F X.obj Y.root))) (pfDiv f) :=
  Pf.nsmul_divBy _ _

variable {P F} in
/-- ★根を**下げる**向きの零因子の変化(掛け算に直した形)。 -/
theorem pfDiv_rtRootIso (A B : C) {dA dB e tA tB : ℕ+} (hA : tA = e * dA) (hB : tB = e * dB)
    (w : HomPf P F (rtObj P F A tA) (rtObj P F B tB)) :
    ((e : ℕ+) : ℕ) • pfDiv ((rtRootIso P F A B hA hB).hom w)
      = Pf.map (Φ.map (P.Base (rtLift P F A hA))) (pfDiv w) := by
  have h := pfDiv_rootIso (rtLift P F A hA) (rtLift_frobType P F A hA)
    (rtLift P F B hB) (rtLift_frobType P F B hB)
    (by rw [rtLift_degFr, rtLift_degFr]) w
  have he : ((e : ℕ+) : ℕ) = ((P.degFr (rtLift P F A hA) : ℕ+) : ℕ) := by
    rw [rtLift_degFr]
  rw [he]
  exact (congrArg (fun t => ((P.degFr (rtLift P F A hA) : ℕ+) : ℕ) • t) h).trans
    (Pf.nsmul_divBy _ _)

variable {P F} in
/-- ★根を**上げる**向きの零因子の変化(掛け算に直した形)。 -/
theorem pfDiv_rtRootIso_inv (A B : C) {dA dB e tA tB : ℕ+} (hA : tA = e * dA)
    (hB : tB = e * dB) (z : HomPf P F (rtObj P F A dA) (rtObj P F B dB)) :
    ((e : ℕ+) : ℕ) • pfDiv z
      = Pf.map (Φ.map (P.Base (rtLift P F A hA)))
        (pfDiv ((rtRootIso P F A B hA hB).inv z)) := by
  have h := pfDiv_rtRootIso A B hA hB ((rtRootIso P F A B hA hB).inv z)
  rw [Iso.inv_hom_id_apply] at h
  exact h

variable (Φ) in
/-- ★`Φ.map` の合成を **`AddMonoidHom` の等式**として。 -/
theorem phi_map_comp_hom {A B E : D} (β : A ⟶ B) (α : B ⟶ E) :
    (Φ.map β).comp (Φ.map α) = Φ.map (β ≫ α) := by
  ext x
  exact (Φ.map_comp α β x).symm

variable (Φ) in
/-- ★`Pf.map ∘ Φ.map` の合成則。 -/
theorem Pf_map_phi_comp {A B E : D} (β : A ⟶ B) (α : B ⟶ E) (x : Pf (Φ.val E)) :
    Pf.map (Φ.map β) (Pf.map (Φ.map α) x) = Pf.map (Φ.map (β ≫ α)) x := by
  rw [Pf.map_map, phi_map_comp_hom]

variable {P F} in
/-- ★★**根を上げたときの `rootBase` の移動**(`rootBase_comp` の中核を独立させたもの)。 -/
theorem rootBase_shift {X Y : PfRootObj P F} (f : HomRoot P F X Y) {tA tB c : ℕ+}
    (hA : tA = c * Y.root) (hB : tB = c * X.root) :
    rootBase f ≫ P.Base (rtExt P F Y.obj tB)
      = P.Base (rtExt P F X.obj tA)
        ≫ pfBase ((rtRootIso P F X.obj Y.obj hA hB).inv f) := by
  rw [← rtLift_ext P F Y.obj hB, P.Base_comp, ← Category.assoc, rootBase_spec,
    Category.assoc, pfBase_rtRootIso_inv, ← Category.assoc, ← P.Base_comp, rtLift_ext]

variable {P F} in
/-- ★★★**根つきでも零因子は合成を保つ**。

★★共通の倍率 `N = X.root · Y.root · Z.root` を掛けると **3 項すべての係数が 1** になり、
`Pf` が perfect なので最後に `N • (−)` の単射性で割る。 -/
theorem rootDiv_comp {X Y Z : PfRootObj P F} (f : HomRoot P F X Y) (g : HomRoot P F Y Z) :
    rootDiv (compRoot P F f g)
      = Pf.map (Φ.map (rootBase f)) (rootDiv g) + ((rootDeg g : ℕ+) : ℕ) • rootDiv f := by
  have hL : ((X.root * Y.root * Z.root : ℕ+) : ℕ) • rootDiv (compRoot P F f g)
      = Pf.map (Φ.map (P.Base (rtExt P F X.obj (Z.root * Y.root))))
        (pfDiv (compPf P F
          ((rtRootIso P F X.obj Y.obj
            (show Z.root * Y.root = Z.root * Y.root from rfl)
            (show Z.root * X.root = Z.root * X.root from rfl)).inv f)
          ((rtRootIso P F Y.obj Z.obj
            (show Z.root * X.root = X.root * Z.root from mul_comm _ _)
            (show Y.root * X.root = X.root * Y.root from mul_comm _ _)).inv g))) := by
    rw [show ((X.root * Y.root * Z.root : ℕ+) : ℕ)
        = ((Y.root : ℕ+) : ℕ) * ((X.root * Z.root : ℕ+) : ℕ) from by push_cast; ring,
      ← smul_smul, rootDiv_spec, ← map_nsmul]
    simp only [compRoot]
    rw [pfDiv_rtRootIso, Pf_map_phi_comp, ← P.Base_comp, rtLift_ext]
  have hT1 : ((X.root * Y.root * Z.root : ℕ+) : ℕ)
        • Pf.map (Φ.map (rootBase f)) (rootDiv g)
      = Pf.map (Φ.map (P.Base (rtExt P F X.obj (Z.root * Y.root))))
        (Pf.map (Φ.map (pfBase ((rtRootIso P F X.obj Y.obj
            (show Z.root * Y.root = Z.root * Y.root from rfl)
            (show Z.root * X.root = Z.root * X.root from rfl)).inv f)))
          (pfDiv ((rtRootIso P F Y.obj Z.obj
            (show Z.root * X.root = X.root * Z.root from mul_comm _ _)
            (show Y.root * X.root = X.root * Y.root from mul_comm _ _)).inv g))) := by
    have hgin : ((X.root * Y.root * Z.root : ℕ+) : ℕ) • rootDiv g
        = Pf.map (Φ.map (P.Base (rtExt P F Y.obj (Z.root * X.root))))
          (pfDiv ((rtRootIso P F Y.obj Z.obj
            (show Z.root * X.root = X.root * Z.root from mul_comm _ _)
            (show Y.root * X.root = X.root * Y.root from mul_comm _ _)).inv g)) := by
      rw [show ((X.root * Y.root * Z.root : ℕ+) : ℕ)
          = ((X.root : ℕ+) : ℕ) * ((Y.root * Z.root : ℕ+) : ℕ) from by push_cast; ring,
        ← smul_smul, rootDiv_spec, ← map_nsmul, pfDiv_rtRootIso_inv, Pf_map_phi_comp,
        ← P.Base_comp, rtLift_ext]
    rw [← map_nsmul, hgin, Pf_map_phi_comp, rootBase_shift, ← Pf_map_phi_comp]
  have hT2 : ((X.root * Y.root * Z.root : ℕ+) : ℕ)
        • (((rootDeg g : ℕ+) : ℕ) • rootDiv f)
      = ((rootDeg g : ℕ+) : ℕ) • Pf.map (Φ.map (P.Base (rtExt P F X.obj (Z.root * Y.root))))
        (pfDiv ((rtRootIso P F X.obj Y.obj
            (show Z.root * Y.root = Z.root * Y.root from rfl)
            (show Z.root * X.root = Z.root * X.root from rfl)).inv f)) := by
    rw [smul_comm,
      show ((X.root * Y.root * Z.root : ℕ+) : ℕ)
        = ((Z.root : ℕ+) : ℕ) * ((X.root * Y.root : ℕ+) : ℕ) from by push_cast; ring,
      ← smul_smul, rootDiv_spec, ← map_nsmul, pfDiv_rtRootIso_inv, Pf_map_phi_comp,
      ← P.Base_comp, rtLift_ext]
  apply (Pf.isPerfectMonoid_pf (X.root * Y.root * Z.root)).1
  show ((X.root * Y.root * Z.root : ℕ+) : ℕ) • rootDiv (compRoot P F f g)
      = ((X.root * Y.root * Z.root : ℕ+) : ℕ)
        • (Pf.map (Φ.map (rootBase f)) (rootDiv g) + ((rootDeg g : ℕ+) : ℕ) • rootDiv f)
  have hdeg : pfDeg ((rtRootIso P F Y.obj Z.obj
      (show Z.root * X.root = X.root * Z.root from mul_comm _ _)
      (show Y.root * X.root = X.root * Y.root from mul_comm _ _)).inv g) = rootDeg g :=
    pfDeg_rtRootIso_inv _ _ _ _ g
  rw [smul_add, hL, hT1, hT2, pfDiv_comp, map_add, map_nsmul, hdeg]

@[simp] theorem Pf.divBy_zero {M : Type w} [AddCommMonoid M] (k : ℕ+) :
    Pf.divBy k (0 : Pf M) = 0 := by
  show Pf.divBy k (Pf.mk 0 1) = Pf.mk 0 1
  rw [Pf.divBy_mk]
  exact Pf.sound 1 (by simp)

variable {P F} in
@[simp] theorem rootDeg_id (X : PfRootObj P F) : rootDeg (idRoot P F X) = 1 := by
  show pfDeg (toHomPf (F := F) (𝟙 (rtObj P F X.obj X.root))) = 1
  rw [pfDeg_toHomPf, P.degFr_id]

variable {P F} in
@[simp] theorem rootBase_id (X : PfRootObj P F) :
    rootBase (idRoot P F X) = 𝟙 (P.toElem.obj X.obj).base := by
  haveI : IsIso (P.Base (rtExt P F X.obj X.root)) := (rtExt_frobType P F X.obj X.root).2
  show P.Base (rtExt P F X.obj X.root)
      ≫ pfBase (toHomPf (F := F) (𝟙 (rtObj P F X.obj X.root)))
      ≫ inv (P.Base (rtExt P F X.obj X.root)) = _
  rw [pfBase_toHomPf, P.Base_id, Category.id_comp, IsIso.hom_inv_id]

variable {P F} in
@[simp] theorem rootDiv_id (X : PfRootObj P F) : rootDiv (idRoot P F X) = 0 := by
  show Pf.divBy (X.root * X.root)
      (Pf.map (Φ.map (P.Base (rtExt P F X.obj X.root)))
        (pfDiv (toHomPf (F := F) (𝟙 (rtObj P F X.obj X.root))))) = 0
  rw [pfDiv_toHomPf, P.Div_id]
  show Pf.divBy _ (Pf.map _ (0 : Pf (Φ.val _))) = 0
  rw [map_zero, Pf.divBy_zero]

/-! ## ★6-k. ★★★**関手 `𝒞^pf ⥤ 𝔽_{Φ^pf}`(原文どおりの `𝒞^pf`)**

原文 (FrdI p.59):
> the left; the lower horizontal arrow is induced by the natural morphism of monoids

★★対象が対 `(A,n)` の `pfRootCategory` からの関手。★これが原文の
`Proposition 3.2, (i)` の右の縦矢印そのものである。 -/

/-- ★★★**[FrdI] Proposition 3.2, (i)** —— `𝒞^pf → 𝔽_{Φ^pf}`(根つき版)。 -/
noncomputable def pfRootToElem :
    PfRootObj P F ⥤ ElemFrobCat (Φ.pfOn (phiSharp P)) where
  obj X := ⟨(P.toElem.obj X.obj).base⟩
  map f := ⟨rootBase f, rootDiv f, rootDeg f⟩
  map_id X := by
    refine ElemFrobCat.Hom.ext ?_ ?_ ?_
    · exact rootBase_id X
    · exact rootDiv_id X
    · exact rootDeg_id X
  map_comp f g := by
    refine ElemFrobCat.Hom.ext ?_ ?_ ?_
    · exact rootBase_comp f g
    · exact rootDiv_comp f g
    · exact rootDeg_comp f g

/-! ## ★6-l. 根つき `𝒞^pf` の pre-Frobenioid 構造

★★`totEpiC` と `connectedC` を根つきの圏で示す。★どちらも
**`compRoot` が 3 つの同型を挟むだけ**であることと、
**`(A,n)` から `(A,1)` へ射がある**ことから、`n = 1` 版に落ちる。 -/

variable {P F} in
/-- ★★**根つきでも全射性**(すべての射が epi)。 -/
theorem pfRoot_epi_cancel {X Y Z : PfRootObj P F} (f : HomRoot P F X Y)
    (g h : HomRoot P F Y Z) (hgh : compRoot P F f g = compRoot P F f h) : g = h := by
  have e := congrArg (fun t : HomRoot P F X Z =>
    (rtRootIso P F X.obj Z.obj (show Z.root * Y.root = Y.root * Z.root from mul_comm _ _)
      (show Y.root * X.root = Y.root * X.root from rfl)).inv t) hgh
  simp only [compRoot, Iso.hom_inv_id_apply] at e
  have e2 := homPf_epi_cancel _ _ _ e
  have e3 := congrArg (fun t => (rtRootIso P F Y.obj Z.obj
    (show Z.root * X.root = X.root * Z.root from mul_comm _ _)
    (show Y.root * X.root = X.root * Y.root from mul_comm _ _)).hom t) e2
  simp only [Iso.inv_hom_id_apply] at e3
  exact e3

/-- ★★**根つき `𝒞^pf` は totally epimorphic**。 -/
theorem pfRoot_totEpi : IsTotallyEpimorphic (PfRootObj P F) :=
  fun _ _ f => ⟨fun g h hgh => pfRoot_epi_cancel f g h hgh⟩

variable {P F} in
/-- ★**`(A,n)` から `(A,1)` へ射がある** —— `A^{(1)} ≅ A → A^{(n)}` を送るだけ。 -/
theorem pfRoot_zigzag_to_one (X : PfRootObj P F) :
    Zigzag X (⟨X.obj, 1⟩ : PfRootObj P F) := by
  haveI := isIso_rtExt_one P F X.obj
  exact Zigzag.of_hom (show HomRoot P F X ⟨X.obj, 1⟩ from
    toHomPf (F := F) (inv (rtExt P F X.obj 1) ≫ rtExt P F X.obj X.root))

/-- ★★**根つき `𝒞^pf` は連結**。 -/
theorem pfRoot_isConnected : IsConnected (PfRootObj P F) := by
  haveI := P.connectedC
  haveI : Nonempty (PfRootObj P F) := ⟨⟨(inferInstance : Nonempty C).some, 1⟩⟩
  refine zigzag_isConnected ?_
  intro X Y
  refine (pfRoot_zigzag_to_one X).trans (Zigzag.trans ?_ (pfRoot_zigzag_to_one Y).symm)
  exact zigzag_obj_of_zigzag (toPfRoot P F) (isPreconnected_zigzag (J := C) X.obj Y.obj)

/-- ★★★**[FrdI] Proposition 3.2, (i)** —— **原文どおりの `𝒞^pf` の pre-Frobenioid 構造**。

原文 (FrdI p.59):
> the left; the lower horizontal arrow is induced by the natural morphism of monoids
-/
noncomputable def pfRootPre :
    PreFrobenioid (PfRootObj P F) (Φ.pfOn (phiSharp P)) where
  toElem := pfRootToElem P F
  divisorial A := Pf.isDivisorial' (P.divisorial A)
  totEpiC := pfRoot_totEpi P F
  totEpiD := P.totEpiD
  connectedC := pfRoot_isConnected P F
  connectedD := P.connectedD

/-! ## ★6-m. 根つきの 1-可換図式 —— `𝒞 → 𝒞^pf` の像での計算則

★`toRootHom f = toHomPf (inv (rtExt A 1) ≫ f ≫ rtExt B 1)` を展開する。
★★**`rtExt A 1` は同型かつ Frobenius 型**(次数 `1`、`Div = 0`)なので、
3 量とも `f` のものに戻る。 -/

variable {P F} in
/-- ★同型の逆射の Frobenius 次数。 -/
theorem degFr_inv_eq_one {A B : C} (u : A ⟶ B) [IsIso u] (h : P.degFr u = 1) :
    P.degFr (inv u) = 1 := by
  have e : P.degFr (inv u ≫ u) = P.degFr u * P.degFr (inv u) := P.degFr_comp (inv u) u
  rw [IsIso.inv_hom_id, P.degFr_id, h, one_mul] at e
  exact e.symm

variable {P F} in
/-- ★同型の底も同型。 -/
theorem isIso_Base_of_isIso {A B : C} (u : A ⟶ B) [IsIso u] : IsIso (P.Base u) :=
  ⟨P.Base (inv u), by rw [← P.Base_comp, IsIso.hom_inv_id, P.Base_id],
    by rw [← P.Base_comp, IsIso.inv_hom_id, P.Base_id]⟩

variable {P F} in
/-- ★同型の逆射の零因子は `0`。 -/
theorem Div_inv_eq_zero {A B : C} (u : A ⟶ B) [IsIso u] (hu : P.Div u = 0) :
    P.Div (inv u) = 0 := by
  have e : P.Div (u ≫ inv u)
      = Φ.map (P.Base u) (P.Div (inv u)) + ((P.degFr (inv u) : ℕ+) : ℕ) • P.Div u :=
    P.Div_comp u (inv u)
  rw [IsIso.hom_inv_id, P.Div_id, hu, smul_zero, add_zero] at e
  exact Φ.map_injective (P.Base u) (by rw [map_zero]; exact e.symm)

variable {P F} in
@[simp] theorem rootDeg_toRootHom {A B : C} (f : A ⟶ B) :
    rootDeg (toRootHom (F := F) f) = P.degFr f := by
  haveI := isIso_rtExt_one P F A
  show pfDeg (toHomPf (F := F) (inv (rtExt P F A 1) ≫ f ≫ rtExt P F B 1)) = _
  rw [pfDeg_toHomPf, P.degFr_comp, P.degFr_comp, rtExt_degFr,
    degFr_inv_eq_one (rtExt P F A 1) (rtExt_degFr P F A 1), one_mul, mul_one]

variable {P F} in
@[simp] theorem rootBase_toRootHom {A B : C} (f : A ⟶ B) :
    rootBase (toRootHom (F := F) f) = P.Base f := by
  haveI := isIso_rtExt_one P F A
  refine (rootBase_uniq (toRootHom (F := F) f) (P.Base f) ?_).symm
  show P.Base f ≫ P.Base (rtExt P F B 1)
      = P.Base (rtExt P F A 1) ≫ pfBase (toHomPf (F := F)
        (inv (rtExt P F A 1) ≫ f ≫ rtExt P F B 1))
  conv_rhs =>
    rw [pfBase_toHomPf, P.Base_comp, P.Base_comp, ← Category.assoc, ← P.Base_comp,
      IsIso.hom_inv_id, P.Base_id, Category.id_comp]

variable {P F} in
@[simp] theorem rootDiv_toRootHom {A B : C} (f : A ⟶ B) :
    rootDiv (toRootHom (F := F) f) = Pf.mk (P.Div f) 1 := by
  haveI := isIso_rtExt_one P F A
  have hA0 : P.Div (rtExt P F A 1) = 0 := (rtExt_frobType P F A 1).1.2
  have hB0 : P.Div (rtExt P F B 1) = 0 := (rtExt_frobType P F B 1).1.2
  have hinv0 : P.Div (inv (rtExt P F A 1)) = 0 := Div_inv_eq_zero _ hA0
  have hx : P.Div (inv (rtExt P F A 1) ≫ f ≫ rtExt P F B 1)
      = Φ.map (P.Base (inv (rtExt P F A 1))) (P.Div f) := by
    rw [P.Div_comp, P.Div_comp, hB0, hinv0]
    simp [rtExt_degFr]
  show Pf.divBy (1 * 1)
      (Pf.map (Φ.map (P.Base (rtExt P F A 1)))
        (pfDiv (toHomPf (F := F) (inv (rtExt P F A 1) ≫ f ≫ rtExt P F B 1)))) = _
  rw [pfDiv_toHomPf, hx, Pf.map_mk, ← Φ.map_comp, ← P.Base_comp, IsIso.hom_inv_id,
    P.Base_id, Φ.map_id, Pf.divBy_mk]
  rfl

/-- ★★★**[FrdI] Proposition 3.2, (i)** —— 射の側の 1-可換性(根つき版)。 -/
theorem pfRootSquare_map {A B : C} (φ : A ⟶ B) :
    (toPfRoot P F ⋙ pfRootToElem P F).map φ = (P.toElem ⋙ elemToPfElem P).map φ := by
  refine ElemFrobCat.Hom.ext ?_ ?_ ?_
  · exact rootBase_toRootHom (F := F) φ
  · exact rootDiv_toRootHom (F := F) φ
  · exact rootDeg_toRootHom (F := F) φ

/-- ★★★**[FrdI] Proposition 3.2, (i)** —— **図式の 1-可換性**(原文どおりの `𝒞^pf`)。

原文 (FrdI p.59):
> the left; the lower horizontal arrow is induced by the natural morphism of monoids

★★これで (i) の 4 辺と 1-可換性がすべて揃った。 -/
noncomputable def pfRootSquare :
    toPfRoot P F ⋙ pfRootToElem P F ≅ P.toElem ⋙ elemToPfElem P :=
  NatIso.ofComponents (fun X => Iso.refl ((toPfRoot P F ⋙ pfRootToElem P F).obj X)) (by
    intro A B φ
    show (toPfRoot P F ⋙ pfRootToElem P F).map φ ≫ 𝟙 _
      = 𝟙 _ ≫ (P.toElem ⋙ elemToPfElem P).map φ
    rw [Category.comp_id, Category.id_comp]
    exact pfRootSquare_map P F φ)

/-! ## ★6-n. (ii) —— 射の型の判定

原文 (FrdI p.59):
> pre-step; base-isomorphism; base-identity endomorphism; isomorphism;

★★まず **`pfRootPre` の 3 量が `rootDeg` / `rootBase` / `rootDiv` そのもの**である
ことを橋渡しする。★これで `MorphismTypes.lean` の語彙が `𝒞^pf` に直接使える。 -/

variable {P F} in
@[simp] theorem pfRootPre_degFr {X Y : PfRootObj P F} (f : X ⟶ Y) :
    (pfRootPre P F).degFr f = rootDeg f := rfl

variable {P F} in
@[simp] theorem pfRootPre_Base {X Y : PfRootObj P F} (f : X ⟶ Y) :
    (pfRootPre P F).Base f = rootBase f := rfl

variable {P F} in
@[simp] theorem pfRootPre_Div {X Y : PfRootObj P F} (f : X ⟶ Y) :
    (pfRootPre P F).Div f = rootDiv f := rfl

/-! ### ★★「点ごと」に決まる 3 つ

★★原文は「cofinal な族が X」と言うが、★**次数・底同型・等長の 3 つは
「どの代表でも X」と同値**である —— ★cofinal 性は要らない。
★これは `repDeg` / `repBase` / `repDiv` が**代表元から直接読める**ことによる。 -/

variable {P F} in
/-- ★代表元の底の**特徴づけ**(`inv` を使わない形)。 -/
theorem repBase_spec {A B : C} (Z : IdxPf P F A B) (φ : Z.right.obj.1 ⟶ Z.right.obj.2) :
    repBase Z φ ≫ P.Base Z.hom.hom.2 = P.Base Z.hom.hom.1 ≫ P.Base φ := by
  haveI : IsIso (P.Base Z.hom.hom.2) := Z.hom.property.2.1.2
  show (P.Base Z.hom.hom.1 ≫ P.Base φ ≫ inv (P.Base Z.hom.hom.2))
      ≫ P.Base Z.hom.hom.2 = _
  rw [Category.assoc, Category.assoc, IsIso.inv_hom_id, Category.comp_id]

/-- ★**四角形の両脇が同型なら、上下の同型性は同値**。 -/
theorem isIso_iff_of_sq {W X Y Z : D} (r : W ⟶ X) (v : X ⟶ Z) (u : W ⟶ Y) (t : Y ⟶ Z)
    [IsIso v] [IsIso u] (h : r ≫ v = u ≫ t) : IsIso r ↔ IsIso t := by
  constructor
  · intro hr
    haveI := hr
    haveI : IsIso (u ≫ t) := by rw [← h]; infer_instance
    exact IsIso.of_isIso_comp_left u t
  · intro ht
    haveI := ht
    haveI : IsIso (r ≫ v) := by rw [h]; infer_instance
    exact IsIso.of_isIso_comp_right r v

variable {P F} in
/-- ★★**次数は代表元の次数そのもの**。 -/
theorem rootDeg_mk {X Y : PfRootObj P F}
    (Z : IdxPf P F (rtObj P F X.obj Y.root) (rtObj P F Y.obj X.root))
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) :
    rootDeg (show HomRoot P F X Y from HomPf.mk Z φ) = P.degFr φ :=
  pfDeg_mk Z φ

variable {P F} in
/-- ★★**底が同型 ⟺ 代表元の底が同型**。 -/
theorem rootBase_isIso_iff {X Y : PfRootObj P F}
    (Z : IdxPf P F (rtObj P F X.obj Y.root) (rtObj P F Y.obj X.root))
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) :
    IsIso (rootBase (show HomRoot P F X Y from HomPf.mk Z φ)) ↔ IsIso (P.Base φ) := by
  haveI hz : IsIso (P.Base Z.hom.hom.2) := Z.hom.property.2.1.2
  haveI ha : IsIso (P.Base Z.hom.hom.1) := Z.hom.property.1.2
  haveI heB : IsIso (P.Base (rtExt P F Y.obj X.root)) := (rtExt_frobType P F Y.obj X.root).2
  haveI heA : IsIso (P.Base (rtExt P F X.obj Y.root)) := (rtExt_frobType P F X.obj Y.root).2
  haveI hv : IsIso (P.Base (rtExt P F Y.obj X.root) ≫ P.Base Z.hom.hom.2) :=
    IsIso.comp_isIso' heB hz
  haveI hu : IsIso (P.Base (rtExt P F X.obj Y.root) ≫ P.Base Z.hom.hom.1) :=
    IsIso.comp_isIso' heA ha
  refine isIso_iff_of_sq _ (P.Base (rtExt P F Y.obj X.root) ≫ P.Base Z.hom.hom.2)
    (P.Base (rtExt P F X.obj Y.root) ≫ P.Base Z.hom.hom.1) _ ?_
  rw [← Category.assoc, rootBase_spec, Category.assoc, pfBase_mk, repBase_spec,
    Category.assoc]
  rfl

/-- ★**sharp な単系では `k • x = 0` から `x = 0`**(`k ≥ 1`)。

★★**在庫の言い換えである**(2026-08-17 の自己修正): 中身は
`isTorsionFreeNaive_of_isSharp`(`MonoidVocabulary.lean`)そのもので、
★**一度これを知らずに 10 行で書き直してしまった**。
★「作る前に**主張の形**で grep する」——`nsmul` で引いて見つからなかったが、
在庫は `IsTorsionFreeNaive`(`n • a = 0 → a = 0`)という**別の名前**で持っていた。 -/
theorem nsmul_eq_zero_of_isSharp {M : Type w} [AddCommMonoid M] (h : IsSharp M)
    {k : ℕ+} {x : M} (hx : ((k : ℕ+) : ℕ) • x = 0) : x = 0 :=
  isTorsionFreeNaive_of_isSharp h x ((k : ℕ+) : ℕ) k.property hx

variable {P F} in
/-- ★★**等長 ⟺ 代表元が等長**。

★★`Pf` の中で `0` になるのは「ある `k` 倍が `0`」だが、
★**`Φ` は divisorial ゆえ sharp**なので、それは `0` に他ならない。 -/
theorem rootDiv_eq_zero_iff {X Y : PfRootObj P F}
    (Z : IdxPf P F (rtObj P F X.obj Y.root) (rtObj P F Y.obj X.root))
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) :
    rootDiv (show HomRoot P F X Y from HomPf.mk Z φ) = 0 ↔ P.Div φ = 0 := by
  show Pf.divBy (X.root * Y.root)
      (Pf.map (Φ.map (P.Base (rtExt P F X.obj Y.root))) (pfDiv (HomPf.mk Z φ))) = 0 ↔ _
  rw [pfDiv_mk]
  show Pf.divBy _ (Pf.map _ (Pf.mk (Φ.map (P.Base Z.hom.hom.1) (P.Div φ)) (repRoot Z)))
      = 0 ↔ _
  rw [Pf.map_mk, Pf.divBy_mk, Pf.mk_eq_zero_iff]
  constructor
  · rintro ⟨k, hk⟩
    rw [← map_nsmul, ← map_nsmul] at hk
    have h1 : Φ.map (P.Base Z.hom.hom.1) (((k : ℕ+) : ℕ) • P.Div φ) = 0 :=
      Φ.map_injective (P.Base (rtExt P F X.obj Y.root)) (hk.trans (map_zero _).symm)
    have h2 : ((k : ℕ+) : ℕ) • P.Div φ = 0 :=
      Φ.map_injective (P.Base Z.hom.hom.1) (h1.trans (map_zero _).symm)
    exact nsmul_eq_zero_of_isSharp (phiSharp P _) h2
  · intro h
    exact ⟨1, by rw [h]; simp⟩

/-! ### ★★ (ii) の 4 主張(`MorphismTypes` の語彙で) -/

variable {P F} in
/-- ★★**[FrdI] Proposition 3.2, (ii)** —— **linear** の判定。 -/
theorem isLinear_mk_iff {X Y : PfRootObj P F}
    (Z : IdxPf P F (rtObj P F X.obj Y.root) (rtObj P F Y.obj X.root))
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) :
    IsLinear (pfRootPre P F) (show HomRoot P F X Y from HomPf.mk Z φ) ↔ IsLinear P φ := by
  show rootDeg (show HomRoot P F X Y from HomPf.mk Z φ) = 1 ↔ _
  rw [rootDeg_mk]
  exact Iff.rfl

variable {P F} in
/-- ★★**[FrdI] Proposition 3.2, (ii)** —— **base-isomorphism** の判定。 -/
theorem isBaseIsomorphism_mk_iff {X Y : PfRootObj P F}
    (Z : IdxPf P F (rtObj P F X.obj Y.root) (rtObj P F Y.obj X.root))
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) :
    IsBaseIsomorphism (pfRootPre P F) (show HomRoot P F X Y from HomPf.mk Z φ)
      ↔ IsBaseIsomorphism P φ :=
  rootBase_isIso_iff Z φ

variable {P F} in
/-- ★★**[FrdI] Proposition 3.2, (ii)** —— **isometry** の判定。 -/
theorem isIsometric_mk_iff {X Y : PfRootObj P F}
    (Z : IdxPf P F (rtObj P F X.obj Y.root) (rtObj P F Y.obj X.root))
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) :
    IsIsometric (pfRootPre P F) (show HomRoot P F X Y from HomPf.mk Z φ)
      ↔ IsIsometric P φ :=
  rootDiv_eq_zero_iff Z φ

variable {P F} in
/-- ★★**[FrdI] Proposition 3.2, (ii)** —— **pre-step** の判定。 -/
theorem isPreStep_mk_iff {X Y : PfRootObj P F}
    (Z : IdxPf P F (rtObj P F X.obj Y.root) (rtObj P F Y.obj X.root))
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) :
    IsPreStep (pfRootPre P F) (show HomRoot P F X Y from HomPf.mk Z φ) ↔ IsPreStep P φ :=
  and_congr (isLinear_mk_iff Z φ) (isBaseIsomorphism_mk_iff Z φ)

variable {P F} in
/-- ★★**[FrdI] Proposition 3.2, (ii)** —— **与えられた Frobenius 次数**の判定。 -/
theorem degFr_mk_iff {X Y : PfRootObj P F}
    (Z : IdxPf P F (rtObj P F X.obj Y.root) (rtObj P F Y.obj X.root))
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) (n : ℕ+) :
    (pfRootPre P F).degFr (show HomRoot P F X Y from HomPf.mk Z φ) = n ↔ P.degFr φ = n := by
  show rootDeg (show HomRoot P F X Y from HomPf.mk Z φ) = n ↔ _
  rw [rootDeg_mk]

/-! ## ★6-j. 段 B′ の残り —— `rootBase_comp` / `rootDiv_comp` の導出(紙の上、2026-08-17)

★`X = (A,n)`, `Y = (B,m)`, `Z = (E,k)`、`eA_d := rtExt A d` と書く。
`compRoot f g = α.hom (compPf (β.inv f) (γ.inv g))` の 3 つの同型は:

| 同型 | 型 | 2 本の持ち上げ |
|---|---|---|
| `β` | `Hom^pf(A^{(km)}, B^{(kn)}) ≅ Hom^pf(A^{(m)}, B^{(n)})` | 次数 `k` |
| `γ` | `Hom^pf(B^{(kn)}, E^{(mn)}) ≅ Hom^pf(B^{(k)}, E^{(m)})` | 次数 `n` |
| `α` | `Hom^pf(A^{(km)}, E^{(mn)}) ≅ Hom^pf(A^{(k)}, E^{(n)})` | 次数 `m` |

★★**要になる既存補題は `rtLift_ext`**: `rtExt A d ≫ rtLift A h = rtExt A t`。

### ★`rootBase_comp` の筋(検算済み)

両辺に `Base(eE_n)` と `Base(rtLift E hα)` を後合成して `inv` を消すと:

- **左辺** → `Base(eA_k ≫ rtLift A hα) ≫ pfBase w = Base(eA_{km}) ≫ pfBase w`
  ここで `w = compPf (β.inv f) (γ.inv g)`、`pfBase w = pfBase (β.inv f) ≫ pfBase (γ.inv g)`
- **右辺** → まず `rootBase g ≫ Base(eE_{mn}) = Base(eB_{kn}) ≫ pfBase (γ.inv g)`
  (`inv(Base(eE_m)) ≫ Base(eE_m) = 𝟙` と `rtLift_ext` と `pfBase_rootIso`)
  次に `rootBase f ≫ Base(eB_{kn}) = Base(eA_{km}) ≫ pfBase (β.inv f)`(同様)

★両者が一致する。★**`Base(eE_n)` は同型ゆえ mono なので消してよい。**

### ★`rootDiv_comp` の筋(同型)

★★**正規化を `X.root * Y.root` に直したので、共通の倍率は
`N = X.root · Y.root · Z.root` で 3 項すべての係数が 1 になる**(紙の上で確認):

- `N • rootDiv (compRoot f g) = Pf.map (Φ.map (Base eA_{km})) (pfDiv w)`
- `N • Pf.map (Φ.map (rootBase f)) (rootDiv g)
   = Pf.map (Φ.map (Base eA_{km})) (Pf.map (Φ.map (pfBase (β.inv f))) (pfDiv (γ.inv g)))`
- `N • (rootDeg g) • rootDiv f
   = (rootDeg g) • Pf.map (Φ.map (Base eA_{km})) (pfDiv (β.inv f))`

★あとは `pfDiv_comp` で `pfDiv w` を 2 項に割り、★**`Pf` が perfect なので
`N • (−)` は単射**で割ればよい。★★**2026-08-17 に実装済み**(`rootDiv_comp`)。

### ★★★ (i) は完成した(2026-08-17)

| 主張 | 実装 |
|---|---|
| 上 `𝒞 → 𝒞^pf` | `toPfRoot`(`Def31Pf.lean`) |
| 左 `𝒞 → 𝔽_Φ` | `P.toElem` |
| 右 `𝒞^pf → 𝔽_{Φ^pf}` | ★**`pfRootToElem`** |
| 下 `𝔽_Φ → 𝔽_{Φ^pf}` | ★**`elemToPfElem`**(`Pf.of` が誘導) |
| **1-可換性** | ★**`pfRootSquare`**(成分は `Iso.refl`) |
| **`𝒞^pf` の pre-Frobenioid 構造** | ★★**`pfRootPre`** |

★★`n = 1` 版(`PfCat`)も並行して残してある(`pfToElem` / `pfSquare` / `pfPre`)——
★**原文の `𝒞^pf` は根つきの方**だが、`n = 1` 版は (ii)(iii) で扱いやすい場面がある。

★★**規律どおり、通っていない `rootBase_comp` / `rootDiv_comp` は置いていない。**
-/

/-! ## ★7. 実装の記録と次の一手

★★**`𝔽_{Φ^pf}` の 3 つの合成則がすべて揃った**(2026-08-17):
`pfDeg_comp` / `pfBase_comp` / `pfDiv_comp`。

## ★7-a. `pfDiv_comp` で効いた 2 つの技(測定)

★★**技 1: 目標の「内側」を書き換えようとしない。**
`Φ.map` は `AddMonoidHom` の coercion を経由するので、目標の内側の項に
`rw` / `simp only` を当てようとすると**そこに項があるのに当たらない**。
★**`Pf.mk` の分子と分母を `congrArg₂ Pf.mk` でまとめて `have` にし、
外側で `rw` する**と通る。実際に通った 3 本:

- `hcomp : repBase Z φ ≫ P.Base b = P.Base a ≫ P.Base φ`
  (`show` で `inv` を露出させ、`simp only [Category.assoc, IsIso.inv_hom_id,
  Category.comp_id]` の後 `rfl`)
- `hY`(`hcomp` を `Φ.map_comp` で両側から挟む。`rw` を使わず
  `Eq.trans` と `congrArg (fun t => Φ.map t (P.Div ψ))` の連結)
- `hPf`(`show Pf.map _ (Pf.mk _ _) = _` → `rw [Pf.map_mk]` → `congrArg₂ Pf.mk hY hmn`)

★★**技 2: 分母共通版 `Pf.mk_add_mk_same`(★0)を別に用意する。**
`Pf.mk_add_mk` は分母が共通のときも分母を掛けるので、そのままだと最後に
`Pf.sound 1` で `(n*n) • (a+b) = n • (n•a + n•b)` を**目標の内側で**示す形になり、
`push_cast` + `simp only [one_mul, smul_add, smul_smul, mul_assoc]` では
**閉じなかった**(`smul_add` が当たらない、と報告される)。
★同じ現象は `MonoidVocabulary.lean` の `Pf` の算術にも出ている
(そちらでも `smul_add` が unused と報告される)。

★★**さらに測定**: `Pf.mk_add_mk_same` を `rw` で当てようとすると
「`Pf.mk ?m ?a + Pf.mk ?m' ?a` が見つからない」と言われる——目標に
**文字通りその形があるのに**である。★**`refine Eq.trans ?_ (… ).symm` で
引数を明示して当てる**と通る。★最後の `rfl` も必要だった
(`rw` の自動 `rfl` は透明度が足りない)。

## ★7-b. ★★**現在地(2026-08-17)** —— `n = 1` 版の (i) は完成した

★★**通ったもの**(`sorry` 0):

| 段 | 内容 | 実装 |
|---|---|---|
| A | 零因子の合成則 | `pfDiv_comp`(★6) |
| B | 関手 `𝒞^pf ⥤ 𝔽_{Φ^pf}` | `pfToElem`(★6-c) |
| C | 1-可換図式 | `pfSquare`(★6-d) |
| D | pre-Frobenioid 構造 | `pfPre`(★6-g) |

★★**段 C は「1-可換」どころか厳密可換だった** —— `pfSquare` の成分は `Iso.refl`。

★★**残るのは `n = 1` からの持ち上げ**である(下の ★7-b')。

## ★7-b'''. ★★**`inv` を主張から追い出す**(2026-08-17 の測定)

★★`repBase_pushIdx` を
`repBase (push V) φ = Base a ≫ repBase V φ ≫ inv (Base b)`
の形で書くと、★**`IsIso (P.Base b)` の instance 合成が本文中で落ちる**
(`haveI` でも `[IsIso (P.Base b)]` の instance 束縛でも同じ)。
さらに `@inv _ _ _ _ (P.Base b) hb.2` と明示すると、今度は
★**`IsIso.hom_inv_id_assoc` などの `rw` が「目標にその形があるのに」当たらない**
(「target expression is not type-correct under the `instances` transparency level」)。

★★**抜け道: `inv` を主張から消す。**
`repBase (push V) φ ≫ P.Base b = P.Base a ≫ repBase V φ`
と書けば主張に instance が要らず、本文でも `haveI` が普通に効く。
★`Base b` は同型なのでこの形でも共役を完全に決めており、**情報は失われない**。
★同じ理由で `pfBase_rootIso` も `≫ P.Base b` の形にしてある。

## ★7-b''. 実装上の測定(`totEpiC` を通したときの記録)

★★**添字圏が「細い」ことが 2 度効く**:
1. `HomPf.eq_iff` が返す 2 本の遷移射 `u v : Z ⟶ W` を `idx_hom_ext` で同一視する
2. `HomPf.eq_iff` の `W` を `idx13` の像へ持ち上げたあと、
   「`v ≫ k ≫ idx13.map t₀`」と「`idx13.map t`」を `idx_hom_ext` で同一視する

★★**要らなかったもの**: 「`f` の 2 つの代表を一致させる」段は**不要**だった。
`exists_rep3 f g` と `exists_rep3 f h` の添字の**上界を取るだけ**で、
`f` は左から、`g` は左から、`h` は右から運べば 3 つが同じ `V` の上に乗る。
★★**最初この段を設計に入れていたが、実装中に落とせると分かった**(150 行の見積もりが
実際には 60 行になった)。

## ★7-b'. ★★**どの圏に構造を入れるのか**(2026-08-17、原文で確認)

原文 (FrdI p.57):
> follows: The objects of Cpf are pairs (A, n), where A

★★**原文の `𝒞^pf` の対象は「対 `(A, n)`」**である(「`(A,n)` は `A` の `n` 乗根と
思え」と原文が言う)。★したがって `Definition 3.1, (iii)` の圏は
**`Def31Pf.lean` の `PfRootObj` / `pfRootCategory`** の方であり、
★**`PfCat`(対象が `C` そのもの)は `n = 1` の部分**にすぎない。

★★**`PfCat ⊆ 𝒞^pf` は圏同値ではない**(見立て): `(A,1) ≅ (A^{(n)}, n)` は成り立つが、
逆に `(A,n)` を `(B,1)` の形にするには `A` が `n` 乗 Frobenius 冪である必要がある。

★★**したがって段 B・C・D は `pfRootCategory` を対象にすべきである。**
★上の `pfDeg` / `pfBase` / `pfDiv` はどれも `HomPf` の上の関数で、
`HomRoot` は `HomPf` の 1 つだから**そのまま使える**が、
★**底と零因子の「行き先」を `(A,n)` 側に付け替える**必要がある——
`A ⟶ rtObj A m` は Frobenius 型ゆえ底同型なので、
★**`repBaseOf` と同じ「`IsIso` を仮引数に割る」技**で共役を取ればよい。

## ★7-c. (i) の残り工程(2026-08-17 に**実測**)

原文 (FrdI p.59):
> the left; the lower horizontal arrow is induced by the natural morphism of monoids

★★**原文の仮定は「`𝒞` が Frobenius-isotropic 型」**である
(`IsOfFrobeniusIsotropicType P`)。★上の `pfDeg` / `pfBase` / `pfDiv` の
井戸定義性には**使っていない**——効いてくるのは (ii) と (iii) である(見立て)。

| 段 | 内容 | 重さ |
|---|---|---|
| A | 分母共通版 `Pf.mk a n + Pf.mk b n = Pf.mk (a+b) n` + `pfDiv_comp` | 小 |
| B | 関手 `𝒞^pf ⥤ ElemFrobCat (Φ^pf)` | 中 |
| C | `𝒞 ⥤ 𝒞^pf` と `𝔽_Φ ⥤ 𝔽_{Φ^pf}` の 1-可換性 | 中 |
| D | `PreFrobenioid (PfCat P F) (Φ.pfOn hsh)` を組む | 中 |

★★**段 D の実測** —— `PreFrobenioid` は **6 フィールド**しかない
(21 フィールドなのは `FrobenioidCore` の方で、そちらは **(iii)** の話):

| # | フィールド | 状態 |
|---|---|---|
| 1 | `toElem : 𝒞^pf ⥤ 𝔽_{Φ^pf}` | ★段 B そのもの |
| 2 | `divisorial : (Φ.pfOn hsh).IsDivisorialOn` | ★**新規 3 本**(下記) |
| 3 | `totEpiC : IsTotallyEpimorphic (PfCat P F)` | ★新規 |
| 4 | `totEpiD : IsTotallyEpimorphic D` | ★★**そのまま `P.totEpiD`**(`𝒟` は不変) |
| 5 | `connectedC : IsConnected (PfCat P F)` | ★新規(下記の移送) |
| 6 | `connectedD : IsConnected D` | ★★**そのまま `P.connectedD`** |

★**#2 の内訳**: `IsDivisorial = IsPreDivisorial ∧ IsSharp` で
`IsPreDivisorial = integral ∧ saturated ∧ of-characteristic-type`。
★★**4 つのうち 2 つは既済**である:

| 条件 | 既済か | 根拠 |
|---|---|---|
| `IsSharp (Pf M)` | ★★**既済** | `Pf.isSharp_pf`(`MonoidVocabulary.lean`) |
| `IsOfCharacteristicType (Pf M)` | ★★**既済** | `isOfCharacteristicType_of_isSharp` に上を渡すだけ |
| `IsIntegralMonoid (Pf M)` | ★新規 | 下記 |
| `IsSaturatedMonoid (Pf M)` | ★新規 | 下記 |

★★**紙の上で確認した道筋(2026-08-17)**:

1. **`IsCancelAdd (Pf M)`**(`M` が簡約的なら)——代表元での計算。
   `IsIntegralMonoid` との往復は `isIntegralMonoid_iff_isCancelAdd` が既にある。
2. ★★**`Gp (Pf M)` は捻れ無し**。`Gp M' = AddLocalization ⊤` なので
   `AddLocalization.induction_on` で `t = toGp u - toGp v` と書け、
   `n • t = 0` から `toGp (n•u) = toGp (n•v)`、
   **integral** で `n•u = n•v`、★**`Pf` が perfect**(`isPerfectMonoid_pf` —— `M` に
   条件が要らない)なので `n•` は単射、よって `u = v`、`t = 0`。
3. **`IsSaturatedMonoid (Pf M)`**: `n • a = toGp x` とする。
   ★`Pf M` は perfect なので `∃ y, n • y = x`。すると `n • (a - toGp y) = 0` で、
   2 より `a = toGp y`。★**つまり「perfect + integral ⟹ saturated」**である。

★★**要点**: ★**完全化そのものが saturated 性を生む** —— 原文が
`Pf` を「always perfect」と述べる(§0)ことが、ここで効く。

★これらは `MonoidVocabulary.lean` に属する汎用補題であり、
★**`Prop32.lean` に置くのは仮置き**——本来の置き場所はそちら。

## ★7-d. (ii) の設計(2026-08-17、原文と語彙を突き合わせ済み)

原文 (FrdI p.59):
> pre-step; base-isomorphism; base-identity endomorphism; isomorphism;

★★**10 の射型はすべて我々の語彙に既にある**(`MorphismTypes.lean`):

| # | 原文 | 我々の述語 |
|---|---|---|
| 1 | morphism of Frobenius type | `IsFrobeniusType` |
| 2 | pre-step | `IsPreStep` |
| 3 | base-isomorphism | `IsBaseIsomorphism` |
| 4 | base-identity endomorphism | `IsBaseIdentity` |
| 5 | isomorphism | `IsIso`(mathlib) |
| 6 | pull-back morphism | `IsPullBack` |
| 7 | isometry | `IsIsometric` |
| 8 | co-angular morphism | `IsCoAngular` |
| 9 | LB-invertible morphism | `IsLBInvertible` |
| 10 | morphism of a given Frobenius degree | `P.degFr φ = n` |

★★**`cofinal collection` の読み**: ある `𝒞^pf` の射を定める `𝒞` の射の族は
`IdxPf P F A B` で添字づけられ、★**この添字圏は filtered である**
(`idxPf_isFiltered`、`Def31Pf.lean`)。★filtered な添字での「cofinal な部分族」は
**「どの添字にも、その先で条件を満たす添字がある」**と同値である。

★★**設計上の要点(見立て)**: 各述語 `X` について
**「`X` は遷移射で保たれる」**(`X φ → X (frobTransport … φ)`)を示せば、
★「cofinal な族が `X`」は**「ある 1 つの代表が `X`」**に潰れる。
★したがって (ii) は **10 個の「遷移射での保存」補題**に分解される見込みである。

★★**(ii) は (i) に依存する** —— `𝒞^pf` の側で `IsPreStep` 等を語るには
`𝒞^pf` の pre-Frobenioid 構造(= 段 D)が要るからである。

## ★7-e. (iii) の設計(2026-08-17)

原文 (FrdI p.59):
> (i), is a Frobenioid of perfect and isotropic type. Moreover, there is a natural

★(iii) は 3 主張ある: **(1) Frobenioid であること**(`FrobenioidCore` の 21 フィールド)、
**(2) perfect かつ isotropic 型**、**(3) `𝒞^pf ≃ (𝒞^pf)^pf`**。

★★**perfect は既済** —— `IsOfPerfectType` は `Φ` が perfect であること、
`isPerfectMonoid_pf` が**無条件で**それを与える。

### ★★ 訂正の記録(2026-08-17)

★★一度**「`𝒞^pf` では Frobenius 型射が同型になる」**と書いたが、★**これは誤り**である。
原文 (ii) の帰納極限は **`(A→A′, B→B′)` の「同じ次数の対」**で添字づけられており、
★**始域と終域を同時に伸ばす**——`Frobenius` 型射を割っているのではない。
実際 `pfDeg` は `ℕ+` のままである(`repDeg Z φ = P.degFr φ`)。
★★**完全化されるのは零因子の側**(`Φ` から `Φ^pf` へ)であって、
Frobenius 次数の側ではない。

### ★★★ isotropic 型が出る正しい筋(2026-08-17、紙の上で確認)

★★**鍵は「cofinal に、添字対象 `A′` 自身が isotropic になる」**ことである。

1. 仮定 `IsOfFrobeniusIsotropicType P` は、各 `A` に **Frobenius 型射
   `A ⟶ Dd`(次数 `m`)で `Dd` が isotropic** なものを与える。
2. ★次数 `m` の倍数 `n` に対し、`A ⟶ A′`(次数 `n`)は **`A ⟶ Dd ⟶ A′` と分解する**
   (`F.frobDegSurj` で `Dd ⟶ Dd′`(次数 `n/m`)を取り、`F.frobDegUniq` で
   `Dd′ ≅ A′`)。
3. ★★**`F.isotropicClosed`**(= `FrobenioidCore` の (vii)(b)、
   「isotropic な対象からの射の終域も isotropic」)より **`A′` は isotropic**。
4. ★次数は割り切れの向きに cofinal なので、★**cofinal に `A′` は isotropic** である。

★そのうえで、`f : A ⟶ B` が `𝒞^pf` の isometric pre-step なら、その代表 `φ : A′ ⟶ B′` は
**`𝒞` の isometric pre-step である**——★**ここに cofinal は要らない**:

| `𝒞^pf` の条件 | `𝒞` の代表での条件 | 理由 |
|---|---|---|
| `pfDeg f = 1` | `P.degFr φ = 1` | ★**定義そのもの**(`repDeg Z φ = P.degFr φ`) |
| `pfBase f` が同型 | `P.Base φ` が同型 | `repBase = Base a ≫ Base φ ≫ inv (Base b)` で、`a` `b` は Frobenius 型ゆえ底同型 |
| `pfDiv f = 0` | `P.Div φ = 0` | ★**`Pf.mk_eq_zero_iff`**(`mk m a = 0 ↔ ∃ k, k • m = 0`)＋ `Φ` が sharp(`k•m = 0` なら `m` は可逆元、よって `0`)＋ `Φ.map_injective` |

★よって 4 で `A′` を isotropic に取れば `φ` は `𝒞` の同型、したがって `f` は
`𝒞^pf` の同型である。∎(見立て)

★★**したがって (iii) の isotropic 性は (ii) を経由しない** —— 依存は
**(i) → (iii)** で足り、(ii) は 21 フィールドのうち **co-angular の 5 つ**にだけ効く。

### ★ 21 フィールドの移送(区分け)

| 群 | フィールド | 手 |
|---|---|---|
| 底 | `baseSurj` / `preStepSpan` / `plBkEquiv` | ★**対象は `𝒞` と同一、`𝒟` は不変**。`toPfCat` で送るだけ(`plBkEquiv` のみ中) |
| Frobenius 次数 | `frobDegSurj` / `frobDegUniq` | 存在は送るだけ。一意性は代表元に降ろす |
| co-angular | `coAngularComp` / `coAngularOfPreStep` / `otriFwd` / `otriBwd` / `otriBase` | ★**(ii) 経由**(cofinal で `𝒞` に降ろす)。★**21 のうち (ii) を要るのはここだけ** |
| 分解 | `arbFactor` / `arbFactorUniq` / `pullBackLB` / `preStepMono` / `preStepFactor` / `preStepFactorUniq` / `preStepFactor'` / `preStepFactorUniq'` | 代表元で分解して送る |
| 忠実性 | `faithfulUpToUnits` | 中 |
| isotropic | `isotropicHullExists` / `isotropicClosed` | ★★**isotropic 型が出れば両方タダ**(hull は恒等射) |

★★**要点**: 21 フィールドのうち **isotropic の 2 つは「isotropic 型」から自動**、
**底の 3 つは `𝒟` が不変なのでほぼ自動**。★実質の山は
**co-angular の 5 つ**(= (ii) 依存)と**分解の 8 つ**である。

★**#5 の手**: `PfCat P F` は**定義から `C` そのもの**(`PfCat _P _F : Type u2 := C`)で、
`toPfCat` は**対象について恒等**である。★したがって `C` の zigzag を
`toPfCat` で送るだけでよい。

★**#3 の手**: `IsTotallyEpimorphic C = ∀ A B (f : A ⟶ B), Epi f`
(**すべての射が epi**)。★`Hom^pf` は `Hom_𝒞` の filtered colimit なので、
`compPf f g = compPf f h` を**共通の段で計算**してから `P.totEpiC` で消す。
★`HomPf.eq_iff` が「ある遷移射を当てた後で一致」としか言わないので、
**そこを詰めるのが山**である(見立て)。
-/

end ABC3.Found.FrdI
