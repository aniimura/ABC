import Mathlib.CategoryTheory.Limits.Types.Filtered

/-!
# `Hom` の帰納極限 —— 汎用の層

★★**[FrdI] は同じ形の帰納極限を 2 度使う**:

原文 (FrdI p.56):
> (ii) Write

原文 (FrdI p.82):
> (Birationalization of a Frobenioid I) For A, B

- **`Hom^pf`**(`Definition 3.1, (ii)`、p.56) ——
  添字は `(A,B)𝒞^{bi-Fr}`(**両側**、`A` と `B` から**出る**同次数の Frobenius 型射の対)、
  遷移写像は `Proposition 1.10, (i)` の四角形。
- **`Hom^birat`**(`Proposition 4.4`、p.82) ——
  添字は `𝒞^{coa-pre}_A`(**片側**、`A` へ**入る** co-angular pre-step)、
  遷移写像は**前合成**。

★★**測定(2026-08-17)**: この 2 つは「添字圏の射のクラスだけが違う」のではない ——
★**(1) 両側 / 片側、(2) コスライス / スライス、(3) 四角形 / 前合成**の 3 点で違う。
★**したがって「射のクラスをパラメータにした 1 つの構成」では両方を落とせない。**

★★**共有できるのはここまで** —— 「**filtered な添字圏上の `Type` 値関手の帰納極限**」。
★合成写像と圏構造は各ケース固有に書く(`pf` は次数の揃え、`birat` は
`Proposition 1.11, (vii)` の引き戻し)。
★★**抽象化しすぎると「何を仮定したか」が読めなくなる**ので、ここで止める。

★★**測定(filtered 性はどこで効くか)**: `mk` / `mk_map` / `exists_rep` / `induction`
の 4 本は **filtered 性を使わない**(任意の小圏の余極限で成り立つ)。
★効くのは `sound` / `eq_iff` / `eq_iff_same` の **3 本だけ**で、
いずれも「**共通の上界へ送って一致**」という有向性そのもの。
-/

namespace ABC3.Found.FrdI

open CategoryTheory CategoryTheory.Limits

universe uJ vJ w

/-! ★★**宇宙について(2026-08-17 の測定)**

★原文は添字圏を「**say, some small skeletal subcategory of** `(A,B)𝒞^{bi-Fr}`」と括弧で断る。
★**これは飾りではない** —— 我々の設定では `C : Type u`, `[Category.{v} C]` なので、
コスライス `(A,B)𝒞^{bi-Fr}` の**対象は `Type (max u v)`** に住み、
一方 `Hom_C(A′,B′)` は `Type v` に住む。
★★**`Type v` の中には、この大きさの添字圏上の余極限は一般に存在しない**
(余極限の台は `Σ j, F.obj j` の商で、`Type (max u v)` に住む)。

★**したがって添字圏の大きさに条件を置かず、目標の宇宙 `w` を自由にし、
`[HasColimit F]` をインスタンス引数で受ける** —— mathlib の
`Types.FilteredColimit` 自身がこの一般性で書かれているので、そのまま乗る。
★応用側(`Hom^pf`)では `ULift` で `Type (max u v)` へ持ち上げて使う。 -/

variable {J : Type uJ} [Category.{vJ} J] (F : J ⥤ Type w) [HasColimit F]

/-- ★★**`Hom` の帰納極限**(汎用の層)——
filtered な添字圏 `J` 上の `Type` 値関手 `F` の帰納極限。 -/
abbrev HomColim : Type w := colimit F

namespace HomColim

/-- ★添字 `j` の元から帰納極限の元を作る。 -/
noncomputable def mk (j : J) (x : F.obj j) : HomColim F := colimit.ι F j x

/-- ★★**共通の上界で一致すれば等しい** —— 帰納極限の要。 -/
theorem sound [IsFiltered J] {i j : J} {x : F.obj i} {y : F.obj j} (k : J) (f : i ⟶ k) (g : j ⟶ k)
    (h : F.map f x = F.map g y) : mk F i x = mk F j y :=
  (Types.FilteredColimit.colimit_eq_iff F).mpr ⟨k, f, g, h⟩

/-- ★★**逆向き** —— 等しければ共通の上界で一致する。

★**`J` が filtered であることがここで効く**(有向性)。 -/
theorem eq_iff [IsFiltered J] {i j : J} {x : F.obj i} {y : F.obj j} :
    mk F i x = mk F j y ↔ ∃ (k : J) (f : i ⟶ k) (g : j ⟶ k), F.map f x = F.map g y :=
  Types.FilteredColimit.colimit_eq_iff F

/-- ★遷移写像で送っても帰納極限では同じ元。 -/
@[simp] theorem mk_map {i j : J} (f : i ⟶ j) (x : F.obj i) :
    mk F j (F.map f x) = mk F i x :=
  by
    show colimit.ι F j (F.map f x) = colimit.ι F i x
    rw [← colimit.w F f]
    rfl

/-- ★**代表元が取れる**。 -/
theorem exists_rep (z : HomColim F) : ∃ (j : J) (x : F.obj j), mk F j x = z :=
  Types.jointly_surjective' z

/-- ★**帰納法の原理** —— 代表元について示せばよい。 -/
theorem induction {P : HomColim F → Prop}
    (h : ∀ (j : J) (x : F.obj j), P (mk F j x)) (z : HomColim F) : P z := by
  obtain ⟨j, x, rfl⟩ := exists_rep F z
  exact h j x

/-- ★**同じ添字での等号判定** —— 有向性で 1 本の遷移射に揃う。 -/
theorem eq_iff_same [IsFiltered J] {i : J} {x y : F.obj i} :
    mk F i x = mk F i y ↔ ∃ (j : J) (f : i ⟶ j), F.map f x = F.map f y :=
  Types.FilteredColimit.isColimit_eq_iff' (colimit.isColimit F) x y

/-! ### ★★普遍性 —— 帰納極限**からの**写像

★★これまで `mk` / `sound` / `eq_iff`(極限**への**射と等号判定)しか置いていなかったが、
**極限からの射**(普遍性)も要る。`HomColim F = colimit F` なので
mathlib の `colimit.desc` がそのまま使える。

★用途: **`𝒞^birat` の普遍性** ——
「co-angular pre-step を同型に送る関手は `𝒞^birat` を経由する」。
これがあると `Corollary 4.10` の `Ψ^birat` や
`Proposition 5.5, (ii)` の梱包を**手で組まずに済む**。 -/

/-- ★降下のための余錐。 -/
noncomputable def descCocone {T : Type w} (f : ∀ j : J, F.obj j → T)
    (hf : ∀ {i j : J} (u : i ⟶ j) (x : F.obj i), f j (F.map u x) = f i x) : Cocone F :=
  Cocone.mk T
    { app := fun j => TypeCat.ofHom (f j)
      naturality := fun _ _ u => by
        ext x
        exact hf u x }

/-- ★★★**帰納極限からの写像** —— 各添字での写像が遷移射と両立すればよい。 -/
noncomputable def desc {T : Type w} (f : ∀ j : J, F.obj j → T)
    (hf : ∀ {i j : J} (u : i ⟶ j) (x : F.obj i), f j (F.map u x) = f i x) :
    HomColim F → T :=
  fun z => colimit.desc F (descCocone F f hf) z

@[simp] theorem desc_mk {T : Type w} (f : ∀ j : J, F.obj j → T)
    (hf : ∀ {i j : J} (u : i ⟶ j) (x : F.obj i), f j (F.map u x) = f i x)
    (j : J) (x : F.obj j) : desc F f hf (mk F j x) = f j x := by
  show colimit.desc F (descCocone F f hf) (colimit.ι F j x) = f j x
  rw [← types_comp_apply (colimit.ι F j) (colimit.desc F (descCocone F f hf)),
    colimit.ι_desc]
  rfl

/-- ★★**極限からの 2 つの写像は代表元で決まる**。 -/
theorem desc_ext {T : Type w} (g h : HomColim F → T)
    (hgh : ∀ (j : J) (x : F.obj j), g (mk F j x) = h (mk F j x)) : g = h :=
  funext (induction F (fun j x => hgh j x))

end HomColim

/-! ## ★★2 変数の写像 —— 「細い」filtered 添字圏の場合

★★**用途**: `Hom^pf` の**合成**。`Hom^pf(A,B) × Hom^pf(B,E) → Hom^pf(A,E)` は
「同じ添字圏の上の 3 つの余極限」の間の 2 変数写像である。

★★**「細い(thin)」= 平行 2 射は等しい**、を仮定する。
★これは我々の添字圏では **`𝒞` が totally epimorphic であること**から出る。
★★**細ければ「共通上界を取れば図式は自動的に可換」**なので、
自然性(＝ well-definedness)の証明が**すべて `thin` 1 本に潰れる。**
-/

namespace HomColim

section Bimap

universe uJ' vJ' w'

variable {J : Type uJ} [Category.{vJ} J] [IsFiltered J]
  {X Y Z : J ⥤ Type w} [HasColimit X] [HasColimit Y] [HasColimit Z]

/-- ★★**2 変数写像を作るための材料**を 1 つの構造体にまとめる。

★★**まとめた理由(実装上の測定)**: `thin` / `m` / `hm` を `variable` に置くと、
Lean 4 の「自動で入る section variable」は **`def` と `theorem` で挙動が違い**、
引数の並びが宣言ごとにずれる。★**構造体 1 個にすれば並びが確定する。** -/
structure BiData (X Y Z : J ⥤ Type w) where
  /-- ★添字圏が**細い**(平行 2 射は等しい)。 -/
  thin : ∀ {i j : J} (f g : i ⟶ j), f = g
  /-- ★各段での 2 変数写像。 -/
  m : ∀ j : J, X.obj j → Y.obj j → Z.obj j
  /-- ★遷移写像と可換(自然性)。 -/
  hm : ∀ {i j : J} (f : i ⟶ j) (x : X.obj i) (y : Y.obj i),
    m j (X.map f x) (Y.map f y) = Z.map f (m i x y)

omit [IsFiltered J] in
/-- ★遷移写像の合成(適用形)。★現行 mathlib の `FunctorToTypes.map_comp_apply` は
非推奨なので、必要な形だけをここで用意する。 -/
theorem map_map (G : J ⥤ Type w) {i j k : J} (f : i ⟶ j) (g : j ⟶ k) (x : G.obj i) :
    G.map g (G.map f x) = G.map (f ≫ g) x := by
  rw [G.map_comp]; rfl

omit [HasColimit X] [HasColimit Y] in
/-- ★★**細さの帰結** —— 2 つの元をどの共通上界へ持ち上げても、同じ元になる。

★★**これが 2 変数写像の well-definedness のすべて**である。 -/
theorem bimap_key (d : BiData X Y Z) {i j k k' : J} (x : X.obj i) (y : Y.obj j)
    (f : i ⟶ k) (g : j ⟶ k) (f' : i ⟶ k') (g' : j ⟶ k') :
    mk Z k (d.m k (X.map f x) (Y.map g y)) = mk Z k' (d.m k' (X.map f' x) (Y.map g' y)) := by
  refine sound Z (IsFiltered.max k k') (IsFiltered.leftToMax k k')
    (IsFiltered.rightToMax k k') ?_
  rw [← d.hm, ← d.hm, map_map, map_map, map_map, map_map,
    d.thin (f ≫ IsFiltered.leftToMax k k') (f' ≫ IsFiltered.rightToMax k k'),
    d.thin (g ≫ IsFiltered.leftToMax k k') (g' ≫ IsFiltered.rightToMax k k')]

/-- ★第 2 引数についての余錐。 -/
noncomputable def bimapCocone (d : BiData X Y Z) (i : J) (x : X.obj i) : Cocone Y :=
  Cocone.mk (HomColim Z)
    { app := fun j => TypeCat.ofHom fun y =>
        mk Z (IsFiltered.max i j)
          (d.m _ (X.map (IsFiltered.leftToMax i j) x) (Y.map (IsFiltered.rightToMax i j) y))
      naturality := fun j j' f => by
        ext y
        show mk Z _ (d.m _ (X.map _ x) (Y.map _ (Y.map f y)))
          = mk Z _ (d.m _ (X.map _ x) (Y.map _ y))
        rw [map_map]
        exact bimap_key d x y _ _ _ _ }

/-- ★2 変数写像の**内側**(第 2 引数についての降下)。 -/
noncomputable def bimapRight (d : BiData X Y Z) (i : J) (x : X.obj i) :
    HomColim Y → HomColim Z :=
  fun y => colimit.desc Y (bimapCocone d i x) y

omit [HasColimit X] in
theorem bimapRight_mk (d : BiData X Y Z) (i j : J) (x : X.obj i) (y : Y.obj j) :
    bimapRight d i x (mk Y j y)
      = mk Z (IsFiltered.max i j)
          (d.m _ (X.map (IsFiltered.leftToMax i j) x) (Y.map (IsFiltered.rightToMax i j) y)) := by
  show colimit.desc Y (bimapCocone d i x) (colimit.ι Y j y) = _
  rw [← types_comp_apply (colimit.ι Y j) (colimit.desc Y (bimapCocone d i x)), colimit.ι_desc]
  rfl

omit [HasColimit X] in
/-- ★★**内側は添字の持ち上げで変わらない**(外側の自然性)。 -/
theorem bimapRight_map (d : BiData X Y Z) {i i' : J} (f : i ⟶ i') (x : X.obj i) :
    bimapRight d i' (X.map f x) = bimapRight d i x := by
  funext w
  refine induction (P := fun w => bimapRight d i' (X.map f x) w = bimapRight d i x w) Y
    (fun j y => ?_) w
  rw [bimapRight_mk, bimapRight_mk, map_map]
  exact bimap_key d x y _ _ _ _

/-- ★第 1 引数についての余錐。 -/
noncomputable def bimapCocone₁ (d : BiData X Y Z) : Cocone X :=
  Cocone.mk (HomColim Y → HomColim Z)
    { app := fun i => TypeCat.ofHom fun x => bimapRight d i x
      naturality := fun i i' f => by
        ext x
        exact bimapRight_map d f x }

/-- ★★**2 変数写像** —— 「細い」filtered 添字圏の上の 3 つの余極限の間。 -/
noncomputable def bimap (d : BiData X Y Z) : HomColim X → HomColim Y → HomColim Z :=
  fun x => colimit.desc X (bimapCocone₁ d) x

theorem bimap_mk (d : BiData X Y Z) (i j : J) (x : X.obj i) (y : Y.obj j) :
    bimap d (mk X i x) (mk Y j y)
      = mk Z (IsFiltered.max i j)
          (d.m _ (X.map (IsFiltered.leftToMax i j) x) (Y.map (IsFiltered.rightToMax i j) y)) := by
  have h : bimap d (mk X i x) = bimapRight d i x := by
    show colimit.desc X (bimapCocone₁ d) (colimit.ι X i x) = _
    rw [← types_comp_apply (colimit.ι X i) (colimit.desc X (bimapCocone₁ d)), colimit.ι_desc]
    rfl
  rw [h, bimapRight_mk]

/-- ★★**同じ添字での計算則** —— 実際に使うのはこの形。 -/
theorem bimap_mk_same (d : BiData X Y Z) (i : J) (x : X.obj i) (y : Y.obj i) :
    bimap d (mk X i x) (mk Y i y) = mk Z i (d.m i x y) := by
  rw [bimap_mk]
  refine (bimap_key d x y _ _ (𝟙 i) (𝟙 i)).trans ?_
  rw [X.map_id, Y.map_id]
  rfl

end Bimap

end HomColim

end ABC3.Found.FrdI
