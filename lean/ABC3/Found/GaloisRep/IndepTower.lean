import ABC3.Found.GaloisRep.TowerLift

/-!
# Galois (G2) 第 74 ブロック —— **★★★★両立する基底の塔**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★ここが `T_l E ≅ ℤ_l²` の要である

各層で `E[l^n] ≅ (ℤ/l^n)²` が言えても、それだけでは逆極限は取れない——
**層をまたいで両立する同型**が要る。★そのために**基底を持ち上げる**:

    (a_n, b_n) が A[l^n] の基底、l·a_{n+1} = a_n, l·b_{n+1} = b_n
      ⟹ (a_{n+1}, b_{n+1}) は A[l^{n+1}] の基底

★★独立性の持ち上げ:`i·a' + j·b' = 0` を `l` 倍すると `i·a + j·b = 0`、
よって `l^n ∣ i, j`。★★★`i = l^n i'` と書くと
`i'·(l^{n−1}a) + j'·(l^{n−1}b) = 0` となり、`l^{n−1}a` は位数 `l` だから `l ∣ i'`。
★★★★合わせて `l^{n+1} ∣ i` ——これが持ち上がりの中身である。

## ★★出発点は `n = 1`

`n = 0` からは持ち上がらない(`A[1] = 0` は情報を持たない)。
★出発点の `A[l]` の基底は第 71 ブロック(`addEquiv_zmod_sq`)から取る。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `IndepPair` | ★位数 `n` の独立な対 |
| `exists_indepPair` | ★`A[n] ≃+ (ℤ/n)²` から基底を取り出す |
| `indepPair_lift` | ★★★**独立性は持ち上がる** |
| `exists_indep_tower` | ★★★★**両立する基底の塔** |
-/

namespace ABC3.Found.GaloisRep

universe u

/-- ★★**位数 `n` の独立な対**——`A[n]` の基底になる 2 元。 -/
structure IndepPair {A : Type u} [AddCommGroup A] (n : ℕ) (a b : A) : Prop where
  ha : n • a = 0
  hb : n • b = 0
  hind : ∀ i j : ℤ, i • a + j • b = 0 → ((n : ℤ) ∣ i ∧ (n : ℤ) ∣ j)

/-- ★`A[n] ≃+ (ℤ/n)²` から独立な対が取れる。 -/
theorem exists_indepPair {A : Type u} [AddCommGroup A] (n : ℕ)
    (h : Nonempty ((ZMod n × ZMod n) ≃+ (nsmulHom A n).ker)) :
    ∃ a b : A, IndepPair n a b := by
  obtain ⟨e⟩ := h
  refine ⟨(e (1, 0) : A), (e (0, 1) : A), ?_, ?_, ?_⟩
  · exact (e (1, 0)).2
  · exact (e (0, 1)).2
  · intro i j hij
    have hsub : i • e (1, 0) + j • e (0, 1) = 0 := by
      ext
      push_cast
      exact hij
    have hz : e ((i • (1, 0) : ZMod n × ZMod n) + j • (0, 1)) = 0 := by
      rw [map_add, map_zsmul, map_zsmul]; exact hsub
    have hz2 : (i • (1, 0) : ZMod n × ZMod n) + j • (0, 1) = 0 := by
      have := congrArg e.symm hz
      simpa using this
    have hi : ((i : ZMod n)) = 0 := by
      have := congrArg Prod.fst hz2
      simpa using this
    have hj : ((j : ZMod n)) = 0 := by
      have := congrArg Prod.snd hz2
      simpa using this
    constructor
    · exact (ZMod.intCast_zmod_eq_zero_iff_dvd i n).1 hi
    · exact (ZMod.intCast_zmod_eq_zero_iff_dvd j n).1 hj

/-- ★★★**独立性は持ち上がる**——`l` 倍で落ちる持ち上げは、また独立な対になる。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to -/
theorem indepPair_lift {A : Type u} [AddCommGroup A] (l m : ℕ) {a b a' b' : A}
    (hab : IndepPair (l ^ (m + 1)) a b)
    (ha' : l • a' = a) (hb' : l • b' = b)
    (hoa : (l ^ (m + 2)) • a' = 0) (hob : (l ^ (m + 2)) • b' = 0) :
    IndepPair (l ^ (m + 2)) a' b' := by
  refine ⟨hoa, hob, ?_⟩
  intro i j hij
  have hla' : ((l : ℤ)) • a' = a := by rw [natCast_zsmul]; exact ha'
  have hlb' : ((l : ℤ)) • b' = b := by rw [natCast_zsmul]; exact hb'
  have h1 : i • a + j • b = 0 := by
    rw [← hla', ← hlb', smul_comm i (l : ℤ) a', smul_comm j (l : ℤ) b', ← smul_add, hij,
      smul_zero]
  obtain ⟨hdi, hdj⟩ := hab.hind i j h1
  rw [show (((l ^ (m + 1) : ℕ)) : ℤ) = (l : ℤ) ^ (m + 1) by push_cast; ring] at hdi hdj
  obtain ⟨i', rfl⟩ := hdi
  obtain ⟨j', rfl⟩ := hdj
  have key : ∀ (c : ℤ) (x : A), ((l : ℤ) ^ (m + 1) * c) • x = ((l : ℤ) ^ m * c) • ((l : ℤ) • x) := by
    intro c x
    rw [← mul_smul]
    congr 1
    ring
  have h2 : ((l : ℤ) ^ m * i') • a + ((l : ℤ) ^ m * j') • b = 0 := by
    rw [← hla', ← hlb', ← key, ← key]
    exact hij
  obtain ⟨hd2, hd3⟩ := hab.hind _ _ h2
  rw [show (((l ^ (m + 1) : ℕ)) : ℤ) = (l : ℤ) ^ (m + 1) by push_cast; ring] at hd2 hd3
  constructor
  · rw [show (((l ^ (m + 2) : ℕ)) : ℤ) = (l : ℤ) ^ (m + 2) by push_cast; ring]
    obtain ⟨t, ht⟩ := hd2
    refine ⟨t, ?_⟩
    have hcancel : (l : ℤ) ^ m * i' = (l : ℤ) ^ m * ((l : ℤ) * t) := by
      rw [ht]; ring
    rcases eq_or_ne ((l : ℤ) ^ m) 0 with h0 | h0
    · have hl0 : (l : ℤ) = 0 := by
        by_contra hne
        exact absurd h0 (pow_ne_zero m hne)
      rw [hl0]
      simp
    · have := mul_left_cancel₀ h0 hcancel
      rw [this]; ring
  · rw [show (((l ^ (m + 2) : ℕ)) : ℤ) = (l : ℤ) ^ (m + 2) by push_cast; ring]
    obtain ⟨t, ht⟩ := hd3
    refine ⟨t, ?_⟩
    have hcancel : (l : ℤ) ^ m * j' = (l : ℤ) ^ m * ((l : ℤ) * t) := by
      rw [ht]; ring
    rcases eq_or_ne ((l : ℤ) ^ m) 0 with h0 | h0
    · have hl0 : (l : ℤ) = 0 := by
        by_contra hne
        exact absurd h0 (pow_ne_zero m hne)
      rw [hl0]
      simp
    · have := mul_left_cancel₀ h0 hcancel
      rw [this]; ring

/-- ★★★★**両立する基底の塔**——`l` 倍で一段下の基底に落ちる基底の列。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★これが `T_l E ≅ ℤ_l²` の骨である。 -/
theorem exists_indep_tower {A : Type u} [AddCommGroup A]
    (hfin : ∀ m : ℕ, 1 ≤ m → Finite (nsmulHom A m).ker)
    (hcard : ∀ m : ℕ, 1 ≤ m → Nat.card (nsmulHom A m).ker = m ^ 2)
    (l : ℕ) (hl : 1 < l) :
    ∃ a b : ℕ → A, (∀ n, IndepPair (l ^ n) (a n) (b n)) ∧
      (∀ n, l • a (n + 1) = a n) ∧ (∀ n, l • b (n + 1) = b n) := by
  haveI := hfin l (by omega)
  have hne : Nonempty ((ZMod l × ZMod l) ≃+ (nsmulHom A l).ker) := by
    refine addEquiv_zmod_sq l (by omega) (fun x => by ext; exact x.2) ?_
    intro d hd
    have hd1 : 1 ≤ d := Nat.pos_of_dvd_of_pos hd (by omega)
    rw [Nat.card_congr (kerSubEquiv l d hd).toEquiv]
    exact hcard d hd1
  obtain ⟨a1, b1, hbase⟩ := exists_indepPair (A := A) l hne
  have step : ∀ (n : ℕ) (p : A × A), ∃ q : A × A,
      IndepPair (l ^ (n + 1)) p.1 p.2 →
        (IndepPair (l ^ (n + 2)) q.1 q.2 ∧ l • q.1 = p.1 ∧ l • q.2 = p.2) := by
    intro n p
    by_cases hp : IndepPair (l ^ (n + 1)) p.1 p.2
    · obtain ⟨x, hx1, hx2⟩ := exists_smul_lift hfin hcard l (n + 1) hl p.1 hp.ha
      obtain ⟨y, hy1, hy2⟩ := exists_smul_lift hfin hcard l (n + 1) hl p.2 hp.hb
      exact ⟨(x, y), fun _ => ⟨indepPair_lift l n hp hx1 hy1 hx2 hy2, hx1, hy1⟩⟩
    · exact ⟨(0, 0), fun h => absurd h hp⟩
  choose F hF using step
  set Q : ℕ → A × A := fun n => Nat.rec (a1, b1) (fun m q => F m q) n with hQdef
  have hQ : ∀ n, IndepPair (l ^ (n + 1)) (Q n).1 (Q n).2 := by
    intro n
    induction n with
    | zero => simpa [hQdef, pow_one] using hbase
    | succ m ih => exact (hF m (Q m) ih).1
  have hQ1 : ∀ m, l • (Q (m + 1)).1 = (Q m).1 := fun m => (hF m (Q m) (hQ m)).2.1
  have hQ2 : ∀ m, l • (Q (m + 1)).2 = (Q m).2 := fun m => (hF m (Q m) (hQ m)).2.2
  refine ⟨fun n => Nat.casesOn n 0 (fun m => (Q m).1),
          fun n => Nat.casesOn n 0 (fun m => (Q m).2), ?_, ?_, ?_⟩
  · intro n
    cases n with
    | zero => exact ⟨by simp, by simp, by intro i j _; exact ⟨one_dvd i, one_dvd j⟩⟩
    | succ m => exact hQ m
  · intro n
    cases n with
    | zero =>
      have := (hQ 0).ha
      rw [pow_one] at this
      simpa using this
    | succ m => exact hQ1 m
  · intro n
    cases n with
    | zero =>
      have := (hQ 0).hb
      rw [pow_one] at this
      simpa using this
    | succ m => exact hQ2 m

/-! ## ★出典の紐付け(`.src`) -/

def exists_indep_tower.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Tate 加群——両立する基底の塔)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
