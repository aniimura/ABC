import ABC3.Found.PGC.ContinuousHomCount
import ABC3.Found.PGC.UnramifiedExtension
import ABC3.Found.PGC.AdjoinFieldClosure

/-!
# `N_n(Gal(K^ur/K)) = n` —— 不分岐 Galois 群の連続指標の個数

[pGC] Proposition 1.2 への経路 C(`ResearchPaper/pgc-goal.md`)のノード F3。

`Found/PGC/InertiaKummer.lean`(第 1033)は上界

  `N_n(Γ_K) ≤ n · gcd(n, q−1)`   (`contHomCard_absGal_le_of_card_ker`)

を、**`#ker(Γ_K の指標 → I_K の指標) ≤ n`** という一点だけを仮定に残した形で
用意している。その一点は `Γ_K/I_K ≅ Gal(K^ur/K)` の連続指標の個数であり、
本ファイルはそれを**等号で**決める:

  `#Hom_cont(Gal(K^ur/K), ℤ/n) = n`.

古典的には `Gal(K^ur/K) ≅ Ẑ` から従うが、**`Ẑ` を作らない**。
mathlib に `Ẑ` は無く(2026-09-05 実測)、射影極限を構成しなくても
「有限次不分岐拡大の塔が `(ℕ, ∣)` で添字づけられる」ことだけで足りる。

## 証明の骨組み

`G := Gal(K^ur/K)`、`A := ℤ/n`(乗法的に書く)とし、`P_N ⊆ G` を
「次数 `N` の不分岐中間体を各点固定する部分群」とする。

1. **`G ↠ ℤ/N` で核が `P_N`**(`exists_unramified_surjective_zmod`)。
   `K^ur` の中の次数 `N` の不分岐中間体 `M_N` を取り、`restrictNormalHom M_N` を
   `Gal(M_N/K) ≅ ℤ/N`(巡回)と合成する。核は `IntermediateField.restrictNormalHom_ker`
   でちょうど `M_N.fixingSubgroup`。
2. **開部分群は必ずどれかの `P_N` を含む**(`exists_unramified_fixingSubgroup_le`)。
   Krull 位相の基本近傍 `E'.fixingSubgroup`(`E'` は `K` 上有限次)を取り、
   `E'` の生成元の**有限集合**を不分岐単項拡大の有向族に押し込む
   (`directed_isUnramifiedAdjoin` + `Directed.finset_le`)。
3. **数え上げ**。`n ∣ N` のとき `#{f : G →* A | P_N ≤ ker f} = #Hom(ℤ/N, ℤ/n) = n`
   (`MonoidHom.liftOfSurjective` + `card_monoidHom_zmod`)。
   連続な `f` は `P_M ≤ ker f` なる `M` を持ち、`N := n·M` に取り替えると
   `{f | P_n ≤ ker f} ⊆ {f | P_N ≤ ker f}` は**個数が等しい有限集合の包含**なので
   等号。すなわち連続指標は必ず `P_n` を殺し、その個数は `n` ちょうど。

★塔の「同じ次数の不分岐拡大は一意」(`adjoin_eq_of_isUnramified`)は使わない。
必要なのは `adjoin_le_of_dvd`(次数が割り切れば塔に載る)だけである。

## 逸脱の記録

* 原典(Mochizuki)はこの主張を局所類体論に投げている。本ファイルは
  相互律も `Ẑ` も経由せず、不分岐拡大の塔だけで示す(経路 C の方針どおり)。
* 「連続」は `Found/PGC/ContinuousHomCount.lean` の流儀、すなわち
  **核が開**であることで定義する(`ZMod n` に位相を入れない)。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC

/-! ## `Hom(ℤ/N, ℤ/n)` の個数(`n ∣ N` のとき `n`)

`Multiplicative (ZMod N)` は `ofAdd 1` で生成されるので、準同型はその行き先で
決まる。逆に `n ∣ N` なら `ZMod.castHom` があるので、任意の行き先が実現できる。 -/

/-- `ofAdd 1` は `Multiplicative (ZMod N)` の生成元。 -/
lemma ofAdd_one_pow_val {N : ℕ} [NeZero N] (a : ZMod N) :
    (Multiplicative.ofAdd (1 : ZMod N)) ^ a.val = Multiplicative.ofAdd a := by
  rw [← ofAdd_nsmul]
  congr 1
  rw [nsmul_eq_mul, mul_one, ZMod.natCast_val, ZMod.cast_id]

/-- `Multiplicative (ZMod N)` からの準同型は `ofAdd 1` の行き先で決まる。 -/
lemma monoidHom_zmod_ext {N : ℕ} [NeZero N] {A : Type*} [Monoid A]
    {f g : Multiplicative (ZMod N) →* A}
    (h : f (Multiplicative.ofAdd (1 : ZMod N)) = g (Multiplicative.ofAdd (1 : ZMod N))) :
    f = g := by
  refine MonoidHom.ext fun x => ?_
  obtain ⟨a, rfl⟩ : ∃ a : ZMod N, Multiplicative.ofAdd a = x := ⟨Multiplicative.toAdd x, rfl⟩
  rw [← ofAdd_one_pow_val a, map_pow, map_pow, h]

/-- `ofAdd 1 ↦ a` で定まる `ℤ/N → ℤ/n`(`n ∣ N` なので `ZMod.castHom` が使える)。 -/
def zmodMonoidHomOf {N n : ℕ} (h : n ∣ N) (a : Multiplicative (ZMod n)) :
    Multiplicative (ZMod N) →* Multiplicative (ZMod n) :=
  AddMonoidHom.toMultiplicative
    ((AddMonoidHom.mulRight (Multiplicative.toAdd a)).comp
      (ZMod.castHom h (ZMod n)).toAddMonoidHom)

@[simp] lemma zmodMonoidHomOf_ofAdd_one {N n : ℕ} (h : n ∣ N) (a : Multiplicative (ZMod n)) :
    zmodMonoidHomOf h a (Multiplicative.ofAdd (1 : ZMod N)) = a := by
  show Multiplicative.ofAdd ((ZMod.castHom h (ZMod n)) 1 * Multiplicative.toAdd a) = a
  rw [map_one, one_mul]
  rfl

/-- `n ∣ N` のとき `Hom(ℤ/N, ℤ/n) ≃ ℤ/n`(`f ↦ f(1)`)。 -/
def zmodMonoidHomEquiv {N n : ℕ} [NeZero N] (h : n ∣ N) :
    (Multiplicative (ZMod N) →* Multiplicative (ZMod n)) ≃ Multiplicative (ZMod n) where
  toFun f := f (Multiplicative.ofAdd 1)
  invFun a := zmodMonoidHomOf h a
  left_inv _ := monoidHom_zmod_ext (by rw [zmodMonoidHomOf_ofAdd_one])
  right_inv a := zmodMonoidHomOf_ofAdd_one h a

/-- **`#Hom(ℤ/N, ℤ/n) = n`**(`n ∣ N`、`N ≠ 0`)。

退化の自己検査:`n ∣ N` を落とすと**偽**(一般には `gcd(N, n)` 個)。
`N = 0` を許すと `ZMod 0 = ℤ` で左辺は無限、`Nat.card = 0 ≠ n` になるので**偽**。 -/
theorem card_monoidHom_zmod {N n : ℕ} (hN : N ≠ 0) (h : n ∣ N) :
    Nat.card (Multiplicative (ZMod N) →* Multiplicative (ZMod n)) = n := by
  haveI : NeZero N := ⟨hN⟩
  haveI : NeZero n := ⟨fun hn => hN (Nat.eq_zero_of_zero_dvd (hn ▸ h))⟩
  rw [Nat.card_congr (zmodMonoidHomEquiv h)]
  simp [Nat.card_eq_fintype_card]

/-- `π : G ↠ ℤ/N` の核を殺す `G →* ℤ/n` は(`n ∣ N` のとき)ちょうど `n` 個。 -/
theorem card_subtype_ker_le {G : Type*} [Group G] {N n : ℕ} (hN : N ≠ 0) (hd : n ∣ N)
    (π : G →* Multiplicative (ZMod N)) (hπ : Function.Surjective π) :
    Nat.card {f : G →* Multiplicative (ZMod n) // π.ker ≤ f.ker} = n := by
  rw [Nat.card_congr (MonoidHom.liftOfSurjective π hπ)]
  exact card_monoidHom_zmod hN hd

variable {p : ℕ} [Fact p.Prime]

/-! ## `K^ur` の中の単項拡大と `K.closure` の中の単項拡大

`K^ur` の元 `y` に対し、中間体は二通りある:`K^ur/K` の中間体
`K⟮y⟯` と、`K.closure/K` の中間体 `K⟮y.val⟯`。前者を `K^ur ⊆ K.closure` で
写すと後者になる(`IntermediateField.adjoin_map`)。この対応で
「不分岐拡大の理論(`K.closure` 側)」を「Krull 位相(`K^ur` 側)」に運ぶ。 -/

/-- `K^ur` の中の `K⟮y⟯` を `K.closure` へ写すと `K⟮y⟯` になる。 -/
lemma map_adjoin_singleton_val (K : PAdicLocalField p) (y : ↥(unramifiedClosure K)) :
    (IntermediateField.adjoin K.carrier ({y} : Set ↥(unramifiedClosure K))).map
        (unramifiedClosure K).val
      = IntermediateField.adjoin K.carrier ({(y : K.closure)} : Set K.closure) := by
  rw [IntermediateField.adjoin_map]
  congr 1
  simp

/-- `K^ur` の中の `K⟮y⟯` と `K.closure` の中の `K⟮y⟯` は `K`-同型。 -/
noncomputable def adjoinUnramifiedAlgEquiv (K : PAdicLocalField p) (y : ↥(unramifiedClosure K)) :
    (IntermediateField.adjoin K.carrier ({y} : Set ↥(unramifiedClosure K)))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({(y : K.closure)} : Set K.closure)) :=
  (IntermediateField.equivMap _ (unramifiedClosure K).val).trans
    (IntermediateField.equivOfEq (map_adjoin_singleton_val K y))

/-- `K.closure` の中で `K⟮z⟯` に入る `K^ur` の元は、`K^ur` の中でも `K⟮z⟯` に入る。 -/
lemma mem_adjoin_of_val_mem (K : PAdicLocalField p) (z s : ↥(unramifiedClosure K))
    (h : (s : K.closure) ∈
      IntermediateField.adjoin K.carrier ({(z : K.closure)} : Set K.closure)) :
    s ∈ IntermediateField.adjoin K.carrier ({z} : Set ↥(unramifiedClosure K)) := by
  rw [← map_adjoin_singleton_val K z] at h
  obtain ⟨w, hw, hwv⟩ := h
  exact (Subtype.ext hwv : w = s) ▸ hw

/-- 塔の包含は `K^ur` の中でも成り立つ。 -/
lemma adjoin_le_adjoin_of_val_le (K : PAdicLocalField p) (x z : ↥(unramifiedClosure K))
    (h : IntermediateField.adjoin K.carrier ({(x : K.closure)} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({(z : K.closure)} : Set K.closure)) :
    IntermediateField.adjoin K.carrier ({x} : Set ↥(unramifiedClosure K))
      ≤ IntermediateField.adjoin K.carrier ({z} : Set ↥(unramifiedClosure K)) := by
  rw [IntermediateField.adjoin_le_iff]
  rintro w hw
  rw [Set.mem_singleton_iff] at hw
  subst hw
  exact mem_adjoin_of_val_mem K z w (h (IntermediateField.mem_adjoin_simple_self _ _))

/-! ## `Gal(K^ur/K) ↠ ℤ/N`、核は次数 `N` の不分岐中間体の固定部分群

`Found/PGC/UnramifiedExtension.lean::exists_surjective_unramifiedClosureGal_to_zmod` は
全射だけを述べていて**核が判らない**(したがって連続性も判らない)。ここでは
中間体を `K^ur` の側で取り直し、`IntermediateField.restrictNormalHom_ker` を当てて
核をちょうど `fixingSubgroup` にする。 -/

/-- **★★★★★★★★`Gal(K^ur/K) ↠ ℤ/N` で核は `K^ur` の中の次数 `N` の
不分岐中間体の固定部分群**。

`K` 上次数 `N` の不分岐拡大 `K_N`(`exists_isCyclic_gal`)を `K^ur` の中間体 `M` として
取り直し、`Gal(M/K) ≅ ℤ/N` と `restrictNormalHom M` を合成する。
`IntermediateField.restrictNormalHom_ker` により核は `M.fixingSubgroup`。

退化の自己検査:`N = 0` では `ZMod 0 = ℤ` で全射になりえない(`Gal(K^ur/K)` の
有限商は `ℤ` を持たない)ので `hN` は落とせない。 -/
theorem exists_unramified_surjective_zmod (K : PAdicLocalField p) {N : ℕ} (hN : N ≠ 0) :
    ∃ (y : ↥(unramifiedClosure K))
      (π : (↥(unramifiedClosure K) ≃ₐ[K.carrier] ↥(unramifiedClosure K)) →*
        Multiplicative (ZMod N)),
      IsUnramifiedAdjoin K (y : K.closure) ∧
      Module.finrank K.carrier
        (IntermediateField.adjoin K.carrier ({(y : K.closure)} : Set K.closure)) = N ∧
      Function.Surjective π ∧
      π.ker = (IntermediateField.adjoin K.carrier
        ({y} : Set ↥(unramifiedClosure K))).fixingSubgroup := by
  obtain ⟨x, hrank, hu, hcyc, hcard⟩ := exists_isCyclic_gal K N hN
  have hxmem : x ∈ unramifiedClosure K :=
    adjoin_le_unramifiedClosure K hu (IntermediateField.mem_adjoin_simple_self _ _)
  set y : ↥(unramifiedClosure K) := ⟨x, hxmem⟩ with hy
  set L := IntermediateField.adjoin K.carrier ({x} : Set K.closure) with hL
  set M := IntermediateField.adjoin K.carrier ({y} : Set ↥(unramifiedClosure K)) with hM
  haveI : Normal K.carrier L := normal_of_isUnramifiedAdjoin K x hu
  haveI : CharZero K.carrier :=
    charZero_of_injective_algebraMap (algebraMap ℚ_[p] K.carrier).injective
  haveI : Algebra.IsSeparable K.carrier L := IntermediateField.isSeparable_tower_bot K.carrier L
  haveI : IsGalois K.carrier L := ⟨⟩
  have eLM : M ≃ₐ[K.carrier] L := adjoinUnramifiedAlgEquiv K y
  haveI : IsGalois K.carrier M := IsGalois.of_algEquiv eLM.symm
  haveI : FiniteDimensional K.carrier M :=
    IntermediateField.finiteDimensional_adjoin (fun z _ => Algebra.IsIntegral.isIntegral z)
  have hcycM : IsCyclic (M ≃ₐ[K.carrier] M) := (MulEquiv.isCyclic (AlgEquiv.autCongr eLM)).mpr hcyc
  have hcardM : Nat.card (M ≃ₐ[K.carrier] M) = N :=
    (Nat.card_congr (AlgEquiv.autCongr eLM).toEquiv).trans hcard
  haveI := hcycM
  have e := (hcardM ▸ (zmodCyclicMulEquiv hcycM)).symm
  haveI : Normal K.carrier ↥(unramifiedClosure K) := normal_unramifiedClosure K
  refine ⟨y, (e : _ ≃* Multiplicative (ZMod N)).toMonoidHom.comp
    (AlgEquiv.restrictNormalHom M), hu, hrank, ?_, ?_⟩
  · exact e.surjective.comp (AlgEquiv.restrictNormalHom_surjective (F := K.carrier)
      (K₁ := M) (E := ↥(unramifiedClosure K)))
  · rw [← IntermediateField.restrictNormalHom_ker M]
    ext σ
    simp [MonoidHom.mem_ker]

/-- **★★★★★★★★`Gal(K^ur/K)` の開部分群は不分岐単項拡大の固定部分群を含む**。

Krull 位相の基本近傍は `K` 上有限次の中間体 `E'` の `fixingSubgroup`
(`krullTopology_mem_nhds_one_iff`)。`E'` は有限生成で、その生成元は有限個の
不分岐単項拡大に入る(`mem_unramifiedClosure_iff`)。族が有向
(`directed_isUnramifiedAdjoin`)なので `Directed.finset_le` で一本にまとまる。

退化の自己検査:`hU`(開)を落とすと**偽**——閉でない部分群や、有限指数でも
開でない部分群には基本近傍が入らない。 -/
theorem exists_unramified_fixingSubgroup_le (K : PAdicLocalField p)
    {U : Subgroup (↥(unramifiedClosure K) ≃ₐ[K.carrier] ↥(unramifiedClosure K))}
    (hU : IsOpen (U : Set (↥(unramifiedClosure K) ≃ₐ[K.carrier] ↥(unramifiedClosure K)))) :
    ∃ y : ↥(unramifiedClosure K), IsUnramifiedAdjoin K (y : K.closure) ∧
      (IntermediateField.adjoin K.carrier
        ({y} : Set ↥(unramifiedClosure K))).fixingSubgroup ≤ U := by
  classical
  have hmem : (U : Set (↥(unramifiedClosure K) ≃ₐ[K.carrier] ↥(unramifiedClosure K)))
      ∈ nhds (1 : (↥(unramifiedClosure K) ≃ₐ[K.carrier] ↥(unramifiedClosure K))) :=
    hU.mem_nhds (one_mem U)
  obtain ⟨M', hfin, hsub⟩ :=
    (krullTopology_mem_nhds_one_iff K.carrier ↥(unramifiedClosure K) _).mp hmem
  haveI := hfin
  obtain ⟨t, ht⟩ := (IntermediateField.essFiniteType_iff.mp inferInstance : M'.FG)
  have hchoice : ∀ s : ↥(unramifiedClosure K), ∃ i : {x : K.closure // IsUnramifiedAdjoin K x},
      (s : K.closure) ∈
        IntermediateField.adjoin K.carrier ({(i : K.closure)} : Set K.closure) := by
    intro s
    obtain ⟨x, hx, hs⟩ := (mem_unramifiedClosure_iff K (s : K.closure)).mp s.2
    exact ⟨⟨x, hx⟩, hs⟩
  choose g hg using hchoice
  haveI := nonempty_isUnramifiedAdjoin K
  obtain ⟨z, hz⟩ := (directed_isUnramifiedAdjoin K).finset_le (t.image g)
  have hzmem : (z : K.closure) ∈ unramifiedClosure K :=
    adjoin_le_unramifiedClosure K z.2 (IntermediateField.mem_adjoin_simple_self _ _)
  refine ⟨⟨(z : K.closure), hzmem⟩, z.2, ?_⟩
  have hle : M' ≤ IntermediateField.adjoin K.carrier
      ({(⟨(z : K.closure), hzmem⟩ : ↥(unramifiedClosure K))} : Set ↥(unramifiedClosure K)) := by
    rw [← ht, IntermediateField.adjoin_le_iff]
    intro s hs
    refine mem_adjoin_of_val_mem K _ s ?_
    exact hz (g s) (Finset.mem_image_of_mem g hs) (hg s)
  intro σ hσ
  exact hsub (IntermediateField.fixingSubgroup_le hle hσ)

/-! ## 本体 -/

/-- **★★★★★★★★★★★★★★★★(F3)`N_n(Gal(K^ur/K)) = n`**。

経路 C(`ResearchPaper/pgc-goal.md`)のノード F3。`Gal(K^ur/K) ≅ Ẑ` 相当の内容を
`Ẑ` を作らずに述べたもの:`Ẑ` の連続指標は `ℤ/n` へちょうど `n` 個ある。

`Found/PGC/InertiaKummer.lean::contHomCard_absGal_le_of_card_ker` が仮定に残していた
一点(`#ker(restrictInertia) ≤ n`)は、`Γ_K/I_K ≅ Gal(K^ur/K)` を経由してこの等式に
帰着する。

**証明**:`P_N` を次数 `N` の不分岐中間体の固定部分群とすると

* `P_N` を殺す指標はちょうど `#Hom(ℤ/N, ℤ/n)` 個で、`n ∣ N` なら `n` 個
  (`card_subtype_ker_le`);
* 連続指標は必ずどれかの `P_M` を殺し(`exists_unramified_fixingSubgroup_le`)、
  `N := n·M` に取り替えると `{f | P_n ≤ ker f} ⊆ {f | P_N ≤ ker f}` が
  **同じ個数の有限集合の包含**になるので等号。

したがって連続指標の全体は「`P_n` を殺すもの」の全体に一致し、個数は `n`。

退化の自己検査。

* `n = 0` では `ZMod 0 = ℤ` で、`Hom_cont(Ẑ, ℤ) = 0` は自明群だから左辺は `1`。
  主張は `1 = 0` となり**偽**。`hn` は落とせない。
* 「連続」(= 核が開)を落とすと**偽**:`Ẑ` から `ℤ/n` への抽象群としての準同型は
  (選択公理のもとで)`n > 1` なら非可算個ある。
* `Gal(K^ur/K)` を `Γ_K` に替えると**偽**——`N_n(Γ_K) = n·gcd(n, q−1)` は
  一般に `n` より大きい(在庫 G1:`n ∣ q−1` なら `n²`)。 -/
theorem contHomCard_unramifiedClosureGal (K : PAdicLocalField p) {n : ℕ} (hn : n ≠ 0) :
    contHomCard ((unramifiedClosure K) ≃ₐ[K.carrier] (unramifiedClosure K)) n = n := by
  obtain ⟨y, π, hyu, hyrank, hπsurj, hπker⟩ := exists_unramified_surjective_zmod K hn
  haveI : FiniteDimensional K.carrier
      (IntermediateField.adjoin K.carrier ({y} : Set ↥(unramifiedClosure K))) :=
    IntermediateField.finiteDimensional_adjoin (fun z _ => Algebra.IsIntegral.isIntegral z)
  have hPopen : IsOpen (((IntermediateField.adjoin K.carrier
      ({y} : Set ↥(unramifiedClosure K))).fixingSubgroup :
      Set (↥(unramifiedClosure K) ≃ₐ[K.carrier] ↥(unramifiedClosure K)))) :=
    IntermediateField.fixingSubgroup_isOpen _
  have hcardP : Nat.card {f : (↥(unramifiedClosure K) ≃ₐ[K.carrier] ↥(unramifiedClosure K)) →*
      Multiplicative (ZMod n) // (IntermediateField.adjoin K.carrier
        ({y} : Set ↥(unramifiedClosure K))).fixingSubgroup ≤ f.ker} = n := by
    rw [← hπker]
    exact card_subtype_ker_le hn dvd_rfl π hπsurj
  have hiff : ∀ f : (↥(unramifiedClosure K) ≃ₐ[K.carrier] ↥(unramifiedClosure K)) →*
      Multiplicative (ZMod n),
      f ∈ contHom (↥(unramifiedClosure K) ≃ₐ[K.carrier] ↥(unramifiedClosure K))
          (Multiplicative (ZMod n))
        ↔ (IntermediateField.adjoin K.carrier
          ({y} : Set ↥(unramifiedClosure K))).fixingSubgroup ≤ f.ker := by
    intro f
    constructor
    · intro hf
      rw [mem_contHom] at hf
      obtain ⟨w, hwu, hwle⟩ := exists_unramified_fixingSubgroup_le K hf
      have hMpos : Module.finrank K.carrier
          (IntermediateField.adjoin K.carrier ({(w : K.closure)} : Set K.closure)) ≠ 0 :=
        Module.finrank_pos.ne'
      obtain ⟨z, πz, hzu, hzrank, hzsurj, hzker⟩ :=
        exists_unramified_surjective_zmod K (N := n * Module.finrank K.carrier
          (IntermediateField.adjoin K.carrier ({(w : K.closure)} : Set K.closure)))
          (mul_ne_zero hn hMpos)
      have h1 : IntermediateField.adjoin K.carrier ({(y : K.closure)} : Set K.closure)
          ≤ IntermediateField.adjoin K.carrier ({(z : K.closure)} : Set K.closure) :=
        adjoin_le_of_dvd K _ _ hyu hzu (by rw [hyrank, hzrank]; exact Dvd.intro _ rfl)
      have h2 : IntermediateField.adjoin K.carrier ({(w : K.closure)} : Set K.closure)
          ≤ IntermediateField.adjoin K.carrier ({(z : K.closure)} : Set K.closure) :=
        adjoin_le_of_dvd K _ _ hwu hzu (by rw [hzrank]; exact Dvd.intro_left _ rfl)
      have hPzP := IntermediateField.fixingSubgroup_le (adjoin_le_adjoin_of_val_le K y z h1)
      have hPzw := IntermediateField.fixingSubgroup_le (adjoin_le_adjoin_of_val_le K w z h2)
      have hcardPz : Nat.card {g : (↥(unramifiedClosure K) ≃ₐ[K.carrier]
          ↥(unramifiedClosure K)) →* Multiplicative (ZMod n) //
          (IntermediateField.adjoin K.carrier
            ({z} : Set ↥(unramifiedClosure K))).fixingSubgroup ≤ g.ker} = n := by
        rw [← hzker]
        exact card_subtype_ker_le (mul_ne_zero hn hMpos) ⟨_, rfl⟩ πz hzsurj
      haveI : Finite {g : (↥(unramifiedClosure K) ≃ₐ[K.carrier]
          ↥(unramifiedClosure K)) →* Multiplicative (ZMod n) //
          (IntermediateField.adjoin K.carrier
            ({z} : Set ↥(unramifiedClosure K))).fixingSubgroup ≤ g.ker} :=
        Nat.finite_of_card_ne_zero (by rw [hcardPz]; exact hn)
      have hbij : Function.Bijective (fun g : {g : (↥(unramifiedClosure K) ≃ₐ[K.carrier]
          ↥(unramifiedClosure K)) →* Multiplicative (ZMod n) //
          (IntermediateField.adjoin K.carrier
            ({y} : Set ↥(unramifiedClosure K))).fixingSubgroup ≤ g.ker} =>
          (⟨g.1, le_trans hPzP g.2⟩ : {g : (↥(unramifiedClosure K) ≃ₐ[K.carrier]
          ↥(unramifiedClosure K)) →* Multiplicative (ZMod n) //
          (IntermediateField.adjoin K.carrier
            ({z} : Set ↥(unramifiedClosure K))).fixingSubgroup ≤ g.ker})) := by
        rw [Nat.bijective_iff_injective_and_card]
        refine ⟨fun a b hab => ?_, by rw [hcardP, hcardPz]⟩
        refine Subtype.ext ?_
        have h' := congrArg Subtype.val hab
        exact h'
      obtain ⟨⟨g, hg⟩, hgeq⟩ := hbij.2 ⟨f, le_trans hPzw hwle⟩
      have hgf : g = f := congrArg Subtype.val hgeq
      exact hgf ▸ hg
    · intro hf
      rw [mem_contHom]
      exact Subgroup.isOpen_mono hf hPopen
  show Nat.card (contHom (↥(unramifiedClosure K) ≃ₐ[K.carrier] ↥(unramifiedClosure K))
    (Multiplicative (ZMod n))) = n
  rw [Nat.card_congr (Equiv.subtypeEquivRight hiff)]
  exact hcardP

end ABC3.Found.PGC
