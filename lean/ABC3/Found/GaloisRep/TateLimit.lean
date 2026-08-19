import ABC3.Found.GaloisRep.PadicLimit

/-!
# Galois (G2) 第 76 ブロック —— **★★★★★★捩れの逆極限は `ℤ_l²`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★群論だけで `T_l ≅ ℤ_l²` が出る

仮定は**捩れの個数だけ**である:

    ∀ m ≥ 1,  A[m] は有限で  #A[m] = m²
      ⟹  lim_n A[l^n]  ≃+  ℤ_l × ℤ_l

★楕円曲線であることは一切使わない——第 65 ブロックが与える個数がすべてである。

## ★★★仕掛け

| 段 | 内容 |
|---|---|
| 両立する基底 `(a_n, b_n)` | 第 74 |
| `Φ(α,β)_n := rep_n(α)·a_n + rep_n(β)·b_n` | ★本ブロック |
| 単射 | ★`0 ≤ rep_n < l^n` かつ `l^n ∣ rep_n` ⟹ `rep_n = 0` |
| 全射 | ★★各層の座標を取り、`ℤ_p` の逆極限表示(第 75)で束ねる |

★★★`rep_n(α) := (toZModPow n α).val` は**加法的でない**が、
`l^n • a_n = 0` なので**法 `l^n` で効く分だけ**で足りる——
これが `zsmul_eq_of_dvd` の役割である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `zsmul_eq_of_dvd` | ★消える元への作用は法 `m` でしか効かない |
| `padicRep` ほか | ★`p` 進整数の整数代表とその合同 |
| `indepPair_bijective` | ★★★★独立な対は座標を与える |
| `limTors` | ★★捩れの逆極限 |
| `addEquiv_limTors` | ★★★★★★**`lim A[l^n] ≃+ ℤ_l²`** |
-/

namespace ABC3.Found.GaloisRep

universe u

/-! ## ★法 `m` でしか効かない作用 -/

/-- ★消える元への作用は、法 `m` でしか効かない。 -/
theorem zsmul_eq_of_dvd {A : Type u} [AddCommGroup A] {x : A} {m : ℕ} (hx : m • x = 0) {i j : ℤ}
    (h : ((m : ℤ)) ∣ i - j) : i • x = j • x := by
  have hz : (i - j) • x = 0 := by
    obtain ⟨t, ht⟩ := h
    rw [ht, mul_comm, mul_smul, show ((m : ℤ)) • x = 0 from by rw [natCast_zsmul]; exact hx,
      smul_zero]
  rw [sub_smul] at hz
  exact sub_eq_zero.mp hz

/-! ## ★`p` 進整数の整数代表 -/

/-- ★`α` の第 `n` 段の整数代表。 -/
noncomputable def padicRep (l : ℕ) [Fact l.Prime] (α : ℤ_[l]) (n : ℕ) : ℤ :=
  ((PadicInt.toZModPow n α).val : ℤ)

theorem padicRep_cast (l : ℕ) [Fact l.Prime] (α : ℤ_[l]) (n : ℕ) :
    ((padicRep l α n : ℤ) : ZMod (l ^ n)) = PadicInt.toZModPow n α := by
  simp [padicRep]

theorem padicRep_nonneg (l : ℕ) [Fact l.Prime] (α : ℤ_[l]) (n : ℕ) : 0 ≤ padicRep l α n :=
  Int.natCast_nonneg _

theorem padicRep_lt (l : ℕ) [Fact l.Prime] (α : ℤ_[l]) (n : ℕ) : padicRep l α n < (l : ℤ) ^ n := by
  have := ZMod.val_lt (PadicInt.toZModPow n α)
  simp only [padicRep]
  exact_mod_cast this

/-- ★段を下げても代表は法 `l^m` で変わらない。 -/
theorem padicRep_dvd_sub (l : ℕ) [Fact l.Prime] (α : ℤ_[l]) {m n : ℕ} (h : m ≤ n) :
    ((l : ℤ) ^ m) ∣ (padicRep l α n - padicRep l α m) := by
  have hz : (((padicRep l α n - padicRep l α m : ℤ)) : ZMod (l ^ m)) = 0 := by
    push_cast
    rw [padicRep_cast]
    have hc : ((padicRep l α n : ℤ) : ZMod (l ^ m))
        = ZMod.castHom (pow_dvd_pow l h) (ZMod (l ^ m)) ((padicRep l α n : ℤ) : ZMod (l ^ n)) := by
      rw [map_intCast]
    rw [hc, padicRep_cast]
    rw [show ZMod.castHom (pow_dvd_pow l h) (ZMod (l ^ m)) (PadicInt.toZModPow n α)
        = PadicInt.toZModPow m α from
      RingHom.congr_fun (PadicInt.zmod_cast_comp_toZModPow m n h) α]
    ring
  have := (ZMod.intCast_zmod_eq_zero_iff_dvd _ (l ^ m)).1 hz
  exact_mod_cast this

/-- ★代表は法 `l^n` で加法的である。 -/
theorem padicRep_dvd_add (l : ℕ) [Fact l.Prime] (α β : ℤ_[l]) (n : ℕ) :
    ((l : ℤ) ^ n) ∣ (padicRep l (α + β) n - (padicRep l α n + padicRep l β n)) := by
  have hz : (((padicRep l (α + β) n - (padicRep l α n + padicRep l β n) : ℤ)) : ZMod (l ^ n))
      = 0 := by
    push_cast
    rw [padicRep_cast, padicRep_cast, padicRep_cast, map_add]
    ring
  have := (ZMod.intCast_zmod_eq_zero_iff_dvd _ (l ^ n)).1 hz
  exact_mod_cast this

/-- ★`ZMod (p^n) → ZMod (p^m)` の落とし込みは、代表元を法 `p^m` で保つ。 -/
theorem castHom_val_dvd (p : ℕ) [Fact p.Prime] {m n : ℕ} (h : m ≤ n) (x : ZMod (p ^ n)) :
    ((p : ℤ) ^ m) ∣
      ((x.val : ℤ) - (((ZMod.castHom (pow_dvd_pow p h) (ZMod (p ^ m))) x).val : ℤ)) := by
  have hz : (((x.val : ℤ) - (((ZMod.castHom (pow_dvd_pow p h) (ZMod (p ^ m))) x).val : ℤ) : ℤ)
      : ZMod (p ^ m)) = 0 := by
    push_cast
    have h1 : ((((ZMod.castHom (pow_dvd_pow p h) (ZMod (p ^ m))) x).val : ℕ) : ZMod (p ^ m))
        = (ZMod.castHom (pow_dvd_pow p h) (ZMod (p ^ m))) x := by
      simp [ZMod.natCast_val, ZMod.cast_id]
    have h2 : ((x.val : ℕ) : ZMod (p ^ m))
        = (ZMod.castHom (pow_dvd_pow p h) (ZMod (p ^ m))) ((x.val : ℕ) : ZMod (p ^ n)) := by
      rw [map_natCast]
    have h3 : ((x.val : ℕ) : ZMod (p ^ n)) = x := by simp [ZMod.natCast_val, ZMod.cast_id]
    rw [h1, h2, h3]
    ring
  have := (ZMod.intCast_zmod_eq_zero_iff_dvd _ (p ^ m)).1 hz
  exact_mod_cast this

/-! ## ★★★★独立な対は座標を与える -/

/-- ★★★★**独立な対は座標を与える**——`(i,j) ↦ i·a + j·b` は `(ℤ/n)² → A` の全単射。 -/
theorem indepPair_bijective {A : Type u} [AddCommGroup A] [Finite A] (n : ℕ) (hn : 1 ≤ n)
    {a b : A} (hab : IndepPair n a b) (hcard : Nat.card A = n ^ 2) :
    Function.Bijective
      (fun p : ZMod n × ZMod n => ((p.1.val : ℤ)) • a + ((p.2.val : ℤ)) • b) := by
  haveI : NeZero n := ⟨by omega⟩
  have hinj : Function.Injective
      (fun p : ZMod n × ZMod n => ((p.1.val : ℤ)) • a + ((p.2.val : ℤ)) • b) := by
    rintro ⟨i1, j1⟩ ⟨i2, j2⟩ hpq
    simp only at hpq
    have h0 : ((i1.val : ℤ) - (i2.val : ℤ)) • a + ((j1.val : ℤ) - (j2.val : ℤ)) • b = 0 := by
      rw [sub_smul, sub_smul]
      calc ((i1.val : ℤ) • a - (i2.val : ℤ) • a) + ((j1.val : ℤ) • b - (j2.val : ℤ) • b)
          = ((i1.val : ℤ) • a + (j1.val : ℤ) • b) - ((i2.val : ℤ) • a + (j2.val : ℤ) • b) := by
            abel
        _ = 0 := sub_eq_zero.mpr hpq
    obtain ⟨hd1, hd2⟩ := hab.hind _ _ h0
    have e1 : i1 = i2 := by
      have hz : (((i1.val : ℤ) - (i2.val : ℤ) : ℤ) : ZMod n) = 0 :=
        (ZMod.intCast_zmod_eq_zero_iff_dvd _ n).2 hd1
      push_cast at hz
      have he := sub_eq_zero.mp hz
      simpa using he
    have e2 : j1 = j2 := by
      have hz : (((j1.val : ℤ) - (j2.val : ℤ) : ℤ) : ZMod n) = 0 :=
        (ZMod.intCast_zmod_eq_zero_iff_dvd _ n).2 hd2
      push_cast at hz
      have he := sub_eq_zero.mp hz
      simpa using he
    rw [e1, e2]
  rw [Nat.bijective_iff_injective_and_card]
  refine ⟨hinj, ?_⟩
  rw [hcard, Nat.card_prod, Nat.card_zmod]
  ring

/-! ## ★★★★★★逆極限 -/

/-- ★★**捩れの逆極限** `lim_n A[l^n]`。 -/
def limTors (A : Type u) [AddCommGroup A] (l : ℕ) :
    AddSubgroup (∀ n : ℕ, (nsmulHom A (l ^ n)).ker) where
  carrier := {f | ∀ n : ℕ, l • ((f (n + 1) : A)) = (f n : A)}
  add_mem' := by
    intro x y hx hy n
    simp only [Set.mem_setOf_eq, Pi.add_apply, AddSubgroup.coe_add, smul_add] at *
    rw [hx n, hy n]
  zero_mem' := by intro n; simp
  neg_mem' := by
    intro x hx n
    simp only [Set.mem_setOf_eq, Pi.neg_apply, AddSubgroup.coe_neg, smul_neg] at *
    rw [hx n]

/-- ★★★★★★**捩れの逆極限は `ℤ_l²` である**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★仮定は「`#A[m] = m²`」だけ——楕円曲線であることは使わない。 -/
theorem addEquiv_limTors {A : Type u} [AddCommGroup A]
    (hfin : ∀ m : ℕ, 1 ≤ m → Finite (nsmulHom A m).ker)
    (hcard : ∀ m : ℕ, 1 ≤ m → Nat.card (nsmulHom A m).ker = m ^ 2)
    (l : ℕ) [Fact l.Prime] :
    Nonempty (limTors A l ≃+ (ℤ_[l] × ℤ_[l])) := by
  have hl : 1 < l := (Fact.out : l.Prime).one_lt
  have hpow1 : ∀ n : ℕ, 1 ≤ l ^ n := fun n => Nat.one_le_pow _ _ (by omega)
  have hcastnat : ∀ n : ℕ, ((l ^ n : ℕ) : ℤ) = (l : ℤ) ^ n := by intro n; push_cast; ring
  obtain ⟨a, b, hind, hac, hbc⟩ := exists_indep_tower hfin hcard l hl
  have hkill : ∀ (n : ℕ) (i : ℤ) (x : A), (l ^ n) • x = 0 → (l ^ n) • (i • x) = 0 := by
    intro n i x hx
    rw [smul_comm, hx, smul_zero]
  have hmem : ∀ (α β : ℤ_[l]) (n : ℕ),
      (padicRep l α n • a n + padicRep l β n • b n) ∈ (nsmulHom A (l ^ n)).ker := by
    intro α β n
    rw [mem_ker_nsmulHom, smul_add, hkill n _ _ (hind n).ha, hkill n _ _ (hind n).hb, add_zero]
  have hdvdsub : ∀ (α : ℤ_[l]) (n : ℕ),
      ((l ^ n : ℕ) : ℤ) ∣ (padicRep l α (n + 1) - padicRep l α n) := by
    intro α n
    rw [hcastnat n]
    exact padicRep_dvd_sub l α (Nat.le_succ n)
  have hdvdadd : ∀ (α β : ℤ_[l]) (n : ℕ),
      ((l ^ n : ℕ) : ℤ) ∣ (padicRep l (α + β) n - (padicRep l α n + padicRep l β n)) := by
    intro α β n
    rw [hcastnat n]
    exact padicRep_dvd_add l α β n
  have hcompat : ∀ (α β : ℤ_[l]) (n : ℕ),
      l • (padicRep l α (n + 1) • a (n + 1) + padicRep l β (n + 1) • b (n + 1))
        = padicRep l α n • a n + padicRep l β n • b n := by
    intro α β n
    rw [smul_add, smul_comm (l : ℕ) (padicRep l α (n + 1)) (a (n + 1)), hac n,
      smul_comm (l : ℕ) (padicRep l β (n + 1)) (b (n + 1)), hbc n,
      zsmul_eq_of_dvd (hind n).ha (hdvdsub α n), zsmul_eq_of_dvd (hind n).hb (hdvdsub β n)]
  set Φ : (ℤ_[l] × ℤ_[l]) →+ limTors A l := AddMonoidHom.mk'
    (fun p => ⟨fun n => ⟨padicRep l p.1 n • a n + padicRep l p.2 n • b n, hmem p.1 p.2 n⟩,
      fun n => hcompat p.1 p.2 n⟩)
    (by
      intro p q
      refine Subtype.ext (funext fun n => Subtype.ext ?_)
      show padicRep l (p.1 + q.1) n • a n + padicRep l (p.2 + q.2) n • b n
        = (padicRep l p.1 n • a n + padicRep l p.2 n • b n)
          + (padicRep l q.1 n • a n + padicRep l q.2 n • b n)
      rw [zsmul_eq_of_dvd (hind n).ha (hdvdadd p.1 q.1 n),
        zsmul_eq_of_dvd (hind n).hb (hdvdadd p.2 q.2 n), add_smul, add_smul]
      abel) with hΦ
  have hΦval : ∀ (p : ℤ_[l] × ℤ_[l]) (n : ℕ),
      (((Φ p).1 n : A)) = padicRep l p.1 n • a n + padicRep l p.2 n • b n := fun _ _ => rfl
  refine ⟨(AddEquiv.ofBijective Φ ⟨?_, ?_⟩).symm⟩
  · rw [injective_iff_map_eq_zero]
    rintro ⟨α, β⟩ hz
    have hz' : ∀ n, padicRep l α n • a n + padicRep l β n • b n = 0 := by
      intro n
      have h1 := congrArg (fun (f : limTors A l) => ((f.1 n : A))) hz
      simpa [hΦval] using h1
    have hzero : ∀ (γ : ℤ_[l]) (n : ℕ),
        ((l ^ n : ℕ) : ℤ) ∣ padicRep l γ n → padicRep l γ n = 0 := by
      intro γ n hd
      have h1 := padicRep_nonneg l γ n
      have h2 := padicRep_lt l γ n
      have h3 : padicRep l γ n % ((l : ℤ) ^ n) = 0 := by
        refine Int.emod_eq_zero_of_dvd ?_
        rw [← hcastnat n]; exact hd
      have h4 : padicRep l γ n % ((l : ℤ) ^ n) = padicRep l γ n := Int.emod_eq_of_lt h1 h2
      omega
    have hα : α = 0 := by
      refine PadicInt.ext_of_toZModPow.1 (fun n => ?_)
      have hd := ((hind n).hind _ _ (hz' n)).1
      have h0 := hzero α n hd
      rw [← padicRep_cast l α n, h0]
      simp
    have hβ : β = 0 := by
      refine PadicInt.ext_of_toZModPow.1 (fun n => ?_)
      have hd := ((hind n).hind _ _ (hz' n)).2
      have h0 := hzero β n hd
      rw [← padicRep_cast l β n, h0]
      simp
    rw [hα, hβ]
    rfl
  · intro f
    have hindsub : ∀ n : ℕ, IndepPair (l ^ n)
        (⟨a n, (hind n).ha⟩ : (nsmulHom A (l ^ n)).ker) ⟨b n, (hind n).hb⟩ := by
      intro n
      refine ⟨by ext; exact (hind n).ha, by ext; exact (hind n).hb, ?_⟩
      intro i j hij
      refine (hind n).hind i j ?_
      have h := congrArg Subtype.val hij
      simpa using h
    have hbij : ∀ n : ℕ, Function.Bijective
        (fun p : ZMod (l ^ n) × ZMod (l ^ n) =>
          (p.1.val : ℤ) • (⟨a n, (hind n).ha⟩ : (nsmulHom A (l ^ n)).ker)
            + (p.2.val : ℤ) • (⟨b n, (hind n).hb⟩ : (nsmulHom A (l ^ n)).ker)) := by
      intro n
      haveI := hfin (l ^ n) (hpow1 n)
      exact indepPair_bijective (l ^ n) (hpow1 n) (hindsub n) (hcard _ (hpow1 n))
    choose g hg using fun n => (hbij n).2 (f.1 n)
    have hgA : ∀ n : ℕ, ((g n).1.val : ℤ) • a n + ((g n).2.val : ℤ) • b n = ((f.1 n : A)) := by
      intro n
      have h := congrArg Subtype.val (hg n)
      simpa using h
    have hcast : ∀ n : ℕ,
        (ZMod.castHom (pow_dvd_pow l (Nat.le_succ n)) (ZMod (l ^ n)) ((g (n + 1)).1),
          ZMod.castHom (pow_dvd_pow l (Nat.le_succ n)) (ZMod (l ^ n)) ((g (n + 1)).2)) = g n := by
      intro n
      refine (hbij n).1 ?_
      rw [hg n]
      refine Subtype.ext ?_
      have hd1 : ((l ^ n : ℕ) : ℤ) ∣ (((g (n + 1)).1.val : ℤ)
          - (((ZMod.castHom (pow_dvd_pow l (Nat.le_succ n)) (ZMod (l ^ n)))
              ((g (n + 1)).1)).val : ℤ)) := by
        rw [hcastnat n]
        exact castHom_val_dvd l (Nat.le_succ n) ((g (n + 1)).1)
      have hd2 : ((l ^ n : ℕ) : ℤ) ∣ (((g (n + 1)).2.val : ℤ)
          - (((ZMod.castHom (pow_dvd_pow l (Nat.le_succ n)) (ZMod (l ^ n)))
              ((g (n + 1)).2)).val : ℤ)) := by
        rw [hcastnat n]
        exact castHom_val_dvd l (Nat.le_succ n) ((g (n + 1)).2)
      show ((ZMod.castHom (pow_dvd_pow l (Nat.le_succ n)) (ZMod (l ^ n)) ((g (n + 1)).1)).val
              : ℤ) • a n
          + ((ZMod.castHom (pow_dvd_pow l (Nat.le_succ n)) (ZMod (l ^ n)) ((g (n + 1)).2)).val
              : ℤ) • b n
        = ((f.1 n : A))
      rw [← zsmul_eq_of_dvd (hind n).ha hd1, ← zsmul_eq_of_dvd (hind n).hb hd2,
        ← f.2 n, ← hgA (n + 1), smul_add,
        smul_comm (l : ℕ) (((g (n + 1)).1.val : ℤ)) (a (n + 1)),
        smul_comm (l : ℕ) (((g (n + 1)).2.val : ℤ)) (b (n + 1)), hac n, hbc n]
    obtain ⟨α, hα⟩ := exists_padicInt_of_compat l (fun n => (g n).1)
      (fun n => congrArg Prod.fst (hcast n))
    obtain ⟨β, hβ⟩ := exists_padicInt_of_compat l (fun n => (g n).2)
      (fun n => congrArg Prod.snd (hcast n))
    refine ⟨(α, β), ?_⟩
    refine Subtype.ext (funext fun n => Subtype.ext ?_)
    rw [hΦval]
    rw [show padicRep l α n = ((g n).1.val : ℤ) by simp only [padicRep]; rw [hα n],
      show padicRep l β n = ((g n).2.val : ℤ) by simp only [padicRep]; rw [hβ n]]
    exact hgA n

/-! ## ★出典の紐付け(`.src`) -/

def addEquiv_limTors.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Tate 加群——捩れの逆極限が Z_l²)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
