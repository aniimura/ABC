import ABC3.Found.PGC.AdjoinIntegers
import ABC3.Found.PGC.LubinTateTowerCompatible

/-!
# `𝒪_K^× ≅ lim_n(𝒪_K/π^n𝒪_K)^×`(`sorry` 無し)

節目(5)(射影極限 `Gal(L_π/K)≅𝒪_K^×`)の最終組み立てに要る部品(i)——
完備な離散付値環の単数群が、その商環の単数群の逆極限に一致する、
という古典的な事実——を、`Mathlib.RingTheory.AdicCompletion`の
環・代数の一般論(`AdicCompletion I R`を新たな環として構成する経路)
を経由せず、**単数群だけを狙う自前の構成**で直接確立する。

## 設計判断: なぜ`AdicCompletion`の一般論を経由しないか

`AdicCompletion I M := {f:Πn,M⧸I^n•⊤ // 両立性}`(定義から`lim_n
M/I^n`そのもの)・`AdicCompletion.of_bijective`(`[IsAdicComplete I M]`
から標準写像`M→AdicCompletion I M`が全単射)は既存だが、これを**環**
として`𝒪_K`と同一視し(`AdicCompletion.ofAlgEquiv`らしき定義が
mathlib-indexに見えるが要追加import・未検証)、そこから単数群への
遺伝を示すのは遠回りと判断した。単数群という**知りたい対象だけ**を
`Πn,(𝒪_K/π^n)^×`の中の両立系として直接切り出し、`IsAdicComplete`が
既に提供する2つの事実(`IsHausdorff`=単射性の源・`IsPrecomplete`=
全射性の源)を単数群の言葉に**直接**翻訳する方が、経由する一般論が
軽い。

## 証明の構造

- **単射性**(`unitReductionHom_injective`): `u`が核に入る
  ⟺`u≡1 (mod π^n)`が全`n`で成り立つ⟺`u-1≡0[SMOD(maximalIdeal)^n•⊤]`
  が全`n`で成り立つ——これに`IsHausdorff.haus`(`IsAdicComplete`の
  半分)を直接適用するだけ。
- **全射性**(`unitReductionHom_surjective`): 与えられた両立系`(v_n)`
  の各成分を(単位性を忘れて)`Ideal.Quotient.mk_surjective`で
  `𝒪_K`の元`f n`へ持ち上げ、`v`自身の両立性(`unitReductionTransition`
  との可換性)を`Ideal.Quotient.factor_mk`で言い換えると`(f n)`が
  `IsPrecomplete.prec`の要求する両立系になっていることが分かる。
  極限`L`を得たら、`n=1`の場合(`L mod π=v_1`が仮定より単位=非零)
  から`L∉maximalIdeal`が従い、局所環の標準事実
  (`IsLocalRing.notMem_maximalIdeal`)で`L`自身が`𝒪_K`の単位になる
  ことを結論する。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical

/-- `n≤m`のときの遷移写像`(𝒪_K/π^n)^×→(𝒪_K/π^m)^×`——`Ideal.Quotient.
factor`(`span{π^n}≤span{π^m}`から誘導される環準同型)を`Units.map`
で持ち上げるだけ。 -/
noncomputable def unitReductionTransition {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    {m n : ℕ} (h : m ≤ n) :
    (𝒪[K.carrier] ⧸ Ideal.span ({π ^ n} : Set (𝒪[K.carrier])))ˣ →*
      (𝒪[K.carrier] ⧸ Ideal.span ({π ^ m} : Set (𝒪[K.carrier])))ˣ :=
  Units.map (Ideal.Quotient.factor (S := Ideal.span ({π ^ n} : Set (𝒪[K.carrier])))
    (T := Ideal.span ({π ^ m} : Set (𝒪[K.carrier])))
    (by rw [Ideal.span_singleton_le_span_singleton]; exact pow_dvd_pow π h))

/-- `Πn,(𝒪_K/π^n)^×`の中の、遷移写像に関して両立する族全体のなす
部分群——`lim_n(𝒪_K/π^n)^×`そのもの。 -/
def CompatibleUnits {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π}) :
    Subgroup (∀ n : ℕ, (𝒪[K.carrier] ⧸ Ideal.span ({π ^ n} : Set (𝒪[K.carrier])))ˣ) where
  carrier := {v | ∀ {m n : ℕ} (h : m ≤ n), unitReductionTransition K hπmax h (v n) = v m}
  one_mem' := by intro m n h; simp
  mul_mem' := by
    intro a b ha hb m n h
    simp only [Set.mem_setOf_eq] at *
    rw [Pi.mul_apply, map_mul, ha h, hb h, Pi.mul_apply]
  inv_mem' := by
    intro a ha m n h
    simp only [Set.mem_setOf_eq] at *
    rw [Pi.inv_apply, map_inv, ha h, Pi.inv_apply]

/-- `𝒪_K^×→(𝒪_K/π^n)^×`(標準的な還元写像)。 -/
noncomputable def unitReductionQuotientMap {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (π : 𝒪[K.carrier]) (n : ℕ) :
    (𝒪[K.carrier])ˣ →* (𝒪[K.carrier] ⧸ Ideal.span ({π ^ n} : Set (𝒪[K.carrier])))ˣ :=
  Units.map (Ideal.Quotient.mk (Ideal.span ({π ^ n} : Set (𝒪[K.carrier]))))

/-- `u∈𝒪_K^×`をその還元の族へ送る——`CompatibleUnits`の元であること
は`Ideal.Quotient.factor_mk`(遷移写像と`mk`の可換性)から直ちに従う。 -/
noncomputable def unitReductionFamily {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (u : (𝒪[K.carrier])ˣ) : CompatibleUnits K hπmax :=
  ⟨fun n => unitReductionQuotientMap K π n u, by
    intro m n h
    apply Units.ext
    show (Ideal.Quotient.factor (S := Ideal.span ({π ^ n} : Set (𝒪[K.carrier])))
        (T := Ideal.span ({π ^ m} : Set (𝒪[K.carrier])))
        (by rw [Ideal.span_singleton_le_span_singleton]; exact pow_dvd_pow π h))
      (Ideal.Quotient.mk (Ideal.span ({π ^ n} : Set (𝒪[K.carrier]))) (u : 𝒪[K.carrier])) =
      Ideal.Quotient.mk (Ideal.span ({π ^ m} : Set (𝒪[K.carrier]))) (u : 𝒪[K.carrier])
    rw [Ideal.Quotient.factor_mk]⟩

/-- `unitReductionFamily`は群準同型。 -/
noncomputable def unitReductionHom {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π}) :
    (𝒪[K.carrier])ˣ →* CompatibleUnits K hπmax where
  toFun := unitReductionFamily K hπmax
  map_one' := by
    apply Subtype.ext
    funext n
    show unitReductionQuotientMap K π n 1 = 1
    simp
  map_mul' a b := by
    apply Subtype.ext
    funext n
    show unitReductionQuotientMap K π n (a * b) = unitReductionQuotientMap K π n a * unitReductionQuotientMap K π n b
    simp

/-- ★★★★★★★★★★**`unitReductionHom`は単射**——`IsHausdorff`
(`IsAdicComplete`の一部、既存インスタンス)から直接。 -/
theorem unitReductionHom_injective {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π}) :
    Function.Injective (unitReductionHom K hπmax) := by
  rw [← MonoidHom.ker_eq_bot_iff, eq_bot_iff]
  intro u hu
  simp only [MonoidHom.mem_ker] at hu
  have hu' : ∀ n : ℕ, unitReductionQuotientMap K π n u = 1 := by
    intro n
    have := congrArg (Subtype.val (p := fun v => v ∈ CompatibleUnits K hπmax)) hu
    exact congrFun this n
  have hval : (u : 𝒪[K.carrier]) - 1 = 0 := by
    apply IsHausdorff.haus (I := IsLocalRing.maximalIdeal (𝒪[K.carrier])) inferInstance
    intro n
    rw [SModEq.zero, hπmax, Ideal.span_singleton_pow]
    have h2 := congrArg Units.val (hu' n)
    unfold unitReductionQuotientMap at h2
    simp only [Units.coe_map, Units.val_one] at h2
    have h3 : Ideal.Quotient.mk (Ideal.span ({π ^ n} : Set (𝒪[K.carrier]))) (u : 𝒪[K.carrier]) =
        Ideal.Quotient.mk (Ideal.span ({π ^ n} : Set (𝒪[K.carrier]))) (1 : 𝒪[K.carrier]) := h2
    have h4 := Ideal.Quotient.eq.mp h3
    simpa using h4
  show u = 1
  apply Units.ext
  rw [Units.val_one]
  exact sub_eq_zero.mp hval

/-- ★★★★★★★★★★★★★★★★**`unitReductionHom`は全射**——`IsPrecomplete`
(`IsAdicComplete`の残り半分)で両立系`(f n)`を実際の極限元`L`へ
持ち上げ、`n=1`の場合(仮定の単位性)から`L`自身が`𝒪_K`の単位で
あることを局所環の標準事実で結論する。 -/
theorem unitReductionHom_surjective {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π}) :
    Function.Surjective (unitReductionHom K hπmax) := by
  intro v
  set f : ℕ → 𝒪[K.carrier] := fun n =>
    Function.surjInv (Ideal.Quotient.mk_surjective (I := Ideal.span ({π ^ n} : Set (𝒪[K.carrier]))))
      ((v.1 n : 𝒪[K.carrier] ⧸ Ideal.span ({π ^ n} : Set (𝒪[K.carrier])))) with hf_def
  have hf_eq : ∀ n, Ideal.Quotient.mk (Ideal.span ({π ^ n} : Set (𝒪[K.carrier]))) (f n) =
      (v.1 n : 𝒪[K.carrier] ⧸ Ideal.span ({π ^ n} : Set (𝒪[K.carrier]))) := by
    intro n
    exact Function.surjInv_eq _ _
  have hcompat : ∀ {m n : ℕ}, m ≤ n → f m ≡ f n [SMOD
      (IsLocalRing.maximalIdeal (𝒪[K.carrier])) ^ m • (⊤ : Submodule (𝒪[K.carrier]) (𝒪[K.carrier]))] := by
    intro m n h
    rw [SModEq.sub_mem, hπmax, Ideal.span_singleton_pow]
    simp only [smul_eq_mul, Ideal.mul_top]
    rw [← Ideal.Quotient.eq]
    have hkey : unitReductionTransition K hπmax h (v.1 n) = v.1 m := v.2 h
    have hval := congrArg Units.val hkey
    unfold unitReductionTransition at hval
    simp only [Units.coe_map] at hval
    rw [hf_eq m, ← Ideal.Quotient.factor_mk (S := Ideal.span ({π ^ n} : Set (𝒪[K.carrier])))
      (T := Ideal.span ({π ^ m} : Set (𝒪[K.carrier])))
      (by rw [Ideal.span_singleton_le_span_singleton]; exact pow_dvd_pow π h), hf_eq n]
    exact hval.symm
  obtain ⟨L, hL⟩ := IsAdicComplete.toIsPrecomplete.prec hcompat
  have hLeq : ∀ n, Ideal.Quotient.mk (Ideal.span ({π ^ n} : Set (𝒪[K.carrier]))) L =
      (v.1 n : 𝒪[K.carrier] ⧸ Ideal.span ({π ^ n} : Set (𝒪[K.carrier]))) := by
    intro n
    have hn := hL n
    rw [SModEq.sub_mem, hπmax, Ideal.span_singleton_pow] at hn
    simp only [smul_eq_mul, Ideal.mul_top] at hn
    rw [← Ideal.Quotient.eq] at hn
    rw [← hn, hf_eq n]
  have hL1unit : IsUnit (Ideal.Quotient.mk (Ideal.span ({π ^ 1} : Set (𝒪[K.carrier]))) L) := by
    rw [hLeq 1]
    exact (v.1 1).isUnit
  have hLnotmax : L ∉ IsLocalRing.maximalIdeal (𝒪[K.carrier]) := by
    rw [hπmax, ← pow_one π]
    intro hcontra
    rw [← Ideal.Quotient.eq_zero_iff_mem] at hcontra
    haveI : Nontrivial (𝒪[K.carrier] ⧸ Ideal.span ({π ^ 1} : Set (𝒪[K.carrier]))) := by
      apply Submodule.Quotient.nontrivial_iff.mpr
      rw [pow_one, ← hπmax]
      exact Ideal.IsPrime.ne_top'
    rw [hcontra] at hL1unit
    exact not_isUnit_zero hL1unit
  have hLunit : IsUnit L := IsLocalRing.notMem_maximalIdeal.mp hLnotmax
  refine ⟨hLunit.unit, ?_⟩
  apply Subtype.ext
  funext n
  show unitReductionQuotientMap K π n hLunit.unit = v.1 n
  apply Units.ext
  show Ideal.Quotient.mk (Ideal.span ({π ^ n} : Set (𝒪[K.carrier]))) L = v.1 n
  exact hLeq n

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`𝒪_K^×≅lim_n(𝒪_K/π^n𝒪_K)^×`**——完備な離散付値環の単数群が、その
商環の単数群の逆極限に一致するという古典的な事実。節目(5)(射影極限
`Gal(L_π/K)≅𝒪_K^×`)の最終組み立てに要る部品(i)。 -/
noncomputable def unitReductionEquiv {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π}) :
    (𝒪[K.carrier])ˣ ≃* CompatibleUnits K hπmax :=
  MulEquiv.ofBijective (unitReductionHom K hπmax)
    ⟨unitReductionHom_injective K hπmax, unitReductionHom_surjective K hπmax⟩

/-- **群論の補助事実**: `S=φ.ker`という等式で得られる第一同型定理の
`≃*`(`h▸(quotientKerEquivOfSurjective φ hsurj) : G⧸S≃*H`)は、`mk`の
上では`φ`そのものとして計算できる——`h`を`subst`すれば`rfl`。

★配管(記録): この事実自体は`subst h; rfl`で一瞬で閉じるが、**大域宣言
として**切り出すことが本質的に重要だった。`have`でこの補題を証明した
まま、直後にそれを`principalUnitsQuotientEquiv`側の具体的な`φ,hsurj,h,u`
に適用しようとすると(`exact key φ hsurj h u`のように4引数すべてを
明示しても)、`h`を型検査する時点で`φ`がまだメタ変数のまま残り、
「型`Eq.{u+1} ?m (MonoidHom.ker ?φ)`だが実際は`Eq.{1} ...`」という
食い違いで失敗する——ローカルな`have`束縛の依存関数適用に特有の
エラボレーション順序の癖(具体例は`tools/lean-idioms.md`)。**大域
定理として**切り出す(`theorem ... := by ...`で環境に積む)と、同じ
4引数の適用が一発で通る。 -/
theorem quotientKerEquivOfSurjective_cast_apply {G H : Type*} [CommGroup G] [Group H]
    (φ : G →* H) (hsurj : Function.Surjective φ) {S : Subgroup G} (h : S = φ.ker) (g : G) :
    (h ▸ (QuotientGroup.quotientKerEquivOfSurjective φ hsurj) : G ⧸ S ≃* H)
      (QuotientGroup.mk g) = φ g := by
  subst h
  rfl

/-- **`principalUnitsQuotientEquiv`の自然性**——`(𝒪_K)^×⧸principalUnits(n)`
から`(𝒪_K/π^n)^×`への第一同型定理由来の同型が、`mk`の上では単なる
還元写像`unitReductionQuotientMap`そのものとして計算できる。節目(5)
の最終組み立てに要る部品(ii)への布石: これで`reciprocityMapLimitCompat`
(`(𝒪_K)^×⧸principalUnits(n)`のレベルで述べられている)を
`(𝒪_K/π^n)^×`/`CompatibleUnits`のレベルの主張へ翻訳できる。 -/
theorem principalUnitsQuotientEquiv_apply_mk {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (n : ℕ) (hn : 1 ≤ n) (u : (𝒪[K.carrier])ˣ) :
    principalUnitsQuotientEquiv K hπmax n hn (QuotientGroup.mk u) =
      unitReductionQuotientMap K π n u := by
  simp only [principalUnitsQuotientEquiv, unitReductionQuotientMap]
  exact quotientKerEquivOfSurjective_cast_apply _ _ (principalUnits_eq_ker K π n) u

/-- **`principalUnitsQuotientEquiv`は遷移写像と可換**——`(𝒪_K)^×⧸
principalUnits(n+1)→(𝒪_K)^×⧸principalUnits(n)`(`QuotientGroup.map`、
`principalUnits_succ_le`由来)を先に取ってから`principalUnitsQuotient
Equiv`で`(𝒪_K/π^n)^×`へ落としても、先に`(𝒪_K/π^(n+1))^×`へ落として
から`unitReductionTransition`で`(𝒪_K/π^n)^×`へ落としても同じ。
`QuotientGroup.induction_on`で`mk u`の場合に還元し、両辺を
`principalUnitsQuotientEquiv_apply_mk`で`unitReductionQuotientMap`
の言葉に翻訳した後は`Ideal.Quotient.factor_mk`(遷移写像と`mk`の
可換性、`unitReductionFamily`の両立性証明と同じ道具)で閉じる。
`reciprocityMapLimitCompat`(`(𝒪_K)^×⧸principalUnits`のレベルでの
n跨ぎ両立性)を`CompatibleUnits`(`(𝒪_K/π^n)^×`のレベル)へ翻訳する
最後の橋。 -/
theorem principalUnitsQuotientEquiv_map_eq {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (n : ℕ) (hn : 1 ≤ n) (x : (𝒪[K.carrier])ˣ ⧸ principalUnits K π (n + 1)) :
    principalUnitsQuotientEquiv K hπmax n hn
        (QuotientGroup.map (principalUnits K π (n + 1)) (principalUnits K π n) (MonoidHom.id _)
          (principalUnits_succ_le K π n) x) =
      unitReductionTransition K hπmax (Nat.le_succ n)
        (principalUnitsQuotientEquiv K hπmax (n + 1) (by omega) x) := by
  induction x using QuotientGroup.induction_on with
  | _ u =>
    rw [QuotientGroup.map_mk, principalUnitsQuotientEquiv_apply_mk, principalUnitsQuotientEquiv_apply_mk,
      MonoidHom.id_apply]
    apply Units.ext
    unfold unitReductionTransition unitReductionQuotientMap
    show Ideal.Quotient.mk (Ideal.span ({π ^ n} : Set (𝒪[K.carrier]))) (u : 𝒪[K.carrier]) =
      Ideal.Quotient.factor
        (S := Ideal.span ({π ^ (n + 1)} : Set (𝒪[K.carrier])))
        (T := Ideal.span ({π ^ n} : Set (𝒪[K.carrier])))
        (by rw [Ideal.span_singleton_le_span_singleton]; exact pow_dvd_pow π (Nat.le_succ n))
        (Ideal.Quotient.mk (Ideal.span ({π ^ (n + 1)} : Set (𝒪[K.carrier]))) (u : 𝒪[K.carrier]))
    rw [Ideal.Quotient.factor_mk]

end ABC3.Found.PGC
